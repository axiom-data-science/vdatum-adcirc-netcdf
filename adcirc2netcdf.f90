      program adcirc2netcdf

!-----------------------------------------------------------------------
!  PURPOSE: transform harmonic analysis files from ADCIRC (fort.53)
!  into netcdf format
!-----------------------------------------------------------------------

      real*4, pointer :: ampl(:,:), phase(:,:), freq(:)
      real*4, pointer :: amp37(:,:), pha37(:,:)
      real*4, pointer :: lon(:), lat(:), depth(:)
      real*4 nodefactor,equilarg,freq37(37),wavespeed37(37)
      integer, pointer :: nm(:,:), match(:)
      integer ncid,ne,np,nfreq,nnodes,n,num,nface,ntype,nvrt
      character*10, pointer :: freqname(:)
      character*80 netcdf_file*80, gridtitle*24, tide_analysis*24
      character freqname37(37)*10
      character globalstr(9)*40
      parameter (pi = 3.1415926535897932)

      data globalstr/'Triangular','2D','ADCIRC',&
     & 'Western North Atlantic Ocean Inversion',&
     & 'ADCIRC refined grid inversion',&
     & 'meton:NASUSER/emyers',&
     & 'NOAA/NOS/OCS/CSDL/MMAP','original','edward.myers@noaa.gov'/
      data freqname37/'M2','S2','N2','K1','M4','O1','M6','MK3','S4',&
     & 'MN4','NU2','S6','MU2','2N2','OO1','LAM2','S1','M1','J1',&
     & 'MM','SSA','SA','MSF','MF','RHO','Q1','T2','R2','2Q1','P1',&
     & '2SM2','M3','L2','2MK3','K2','M8','MS4'/
      data wavespeed37/28.9841042,30.0000000,28.4397295,15.0410686,&
     & 57.9682084,13.9430356,86.9523127,44.0251729,60.0000000,&
     & 57.4238337,28.5125831,90.0000000,27.9682084,27.8953548,&
     & 16.1391017,29.4556253,15.0000000,14.4966939,15.5854433,&
     & 0.5443747,0.0821373,0.0410686,1.0158958,1.0980331,13.4715145,&
     & 13.3986609,29.9589333,30.0410667,12.8542862,14.9589314,&
     & 31.0158958,43.4761563,29.5284789,42.9271398,30.0821373,&
     & 115.9364166,58.9841042/

      do i=1,37
        freq37(i)=(wavespeed37(i)*pi/180.0)/3600.0
        freqname37(i)(10:10)=char(0)
      end do

      if (globalstr(1) == "Triangular") then
          nface=3
      else
          nface=4
      endif

!  -  read ADCIRC fort.14 file, in longitude/latitude format
      open(unit=14,file='fort.14',status='old')
      read(14,'(a24)') gridtitle
      read(14,*) ne,np
      allocate(lon(np),lat(np),depth(np))
      allocate(nm(nface,ne))
      do i=1,np
        read(14,*) num,lon(i),lat(i),depth(i)
      end do
      do i=1,ne
        read(14,*) num,ntype,nm(1,i),nm(2,i),nm(3,i)
        if (ntype.ne.nface) then
          write(*,*) 'ntype does not equal nface'
          stop
        endif
      end do
      read(14,*) nope
      read(14,*) neta
      close(14)
      write(*,*) 'Done reading grid'
!  -  read ADCIRC fort.53 file
      open(unit=53,file='fort.53',status='old')
      read(53,*)nfreq
      write(*,*) 'nfreq =',nfreq
      allocate(freqname(nfreq),freq(nfreq))
      allocate(match(nfreq))
      write(*,*) 'np =',np
      allocate(ampl(nfreq,np),phase(nfreq,np))
      allocate(amp37(37,np),pha37(37,np))
3679  format(1x,e20.10,1x,f10.7,1x,f12.8,2x,a10)
      do i=1,nfreq
        read(53,3679)freq(i),nodefactor,equilarg,freqname(i)
        write(*,*) i,freq(i),freqname(i)
      end do
      call match_to_nos37(nfreq,freqname,match)
      read(53,*)nnodes
      write(*,*) 'nnodes=',nnodes
      if (nnodes/=np) then
        write(*,*) 'nnodes does not equal np'
        stop
      endif
      do i=1,nnodes
        do j=1,37
          amp37(j,i)=-99999.0
          pha37(j,i)=-99999.0
        end do
        read(53,*) num
        do j=1,nfreq
          read(53,*) ampl(j,i),phase(j,i)
          if (match(j)/=-1) then
            amp37(match(j),i)=ampl(j,i)
            pha37(match(j),i)=phase(j,i)
          endif
        end do
      end do
      close(53)

! Write tidal information to the netcdf file
      netcdf_file="adcirc53.nc"
      tide_analysis="ADCIRC least squares"
      call get_timeanalysis(neta,time_analysis)
      write(*,*) 'Time of harmonic analysis is ',time_analysis
      nvrt=1
      call write_tidescdf_fem(netcdf_file,ncid,globalstr,&
     & ne,nnodes,nvrt,37,nm,lon,lat,depth,freqname37,freq37,&
     & amp37,pha37,-1.,-1.,-1.,-1.,tide_analysis,time_analysis)

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine match_to_nos37(nn,iname,imatch)

      character*10 iname(nn)
      character frname37(37)*10
      integer imatch(nn)
      data frname37/'M2','S2','N2','K1','M4','O1','M6','MK3','S4',&
     & 'MN4','NU2','S6','MU2','2N2','OO1','LAM2','S1','M1','J1',&
     & 'MM','SSA','SA','MSF','MF','RHO','Q1','T2','R2','2Q1','P1',&
     & '2SM2','M3','L2','2MK3','K2','M8','MS4'/

      do i=1,nn
        imatch(i)=-1
        do j=1,37
          if (iname(i)==frname37(j)) then
            imatch(i)=j
          endif
        end do
      end do

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_timeanalysis(NETA,time_harmonic)

      real*4 time_harmonic


      OPEN(15,FILE='fort.15',status='old')
      READ(15,'(A32)') RUNDES
      READ(15,'(A24)') RUNID
      READ(15,*) NFOVER
      READ(15,*) NABOUT
      READ(15,*) NSCREEN
      READ(15,*) IHOT
      READ(15,*) ICS
      READ(15,*) IM
      READ(15,*) NOLIBF
      READ(15,*) NOLIFA
      READ(15,*) NOLICA
      READ(15,*) NOLICAT
      READ(15,*) NWP
      READ(15,*) NCOR
      READ(15,*) NTIP
      READ(15,*) NWS
      NRS=0
      IF(ABS(NWS).GE.100) THEN
        NRS=1
        NWS=(ABS(NWS)-100)*(NWS/ABS(NWS))
      ENDIF
      READ(15,*) NRAMP
      READ(15,*) G
      READ(15,*) TAU0
      READ(15,*) DTDP
      READ(15,*) STATIM
      READ(15,*) REFTIM
      IF((NWS.EQ.0).AND.(NRS.EQ.1)) READ(15,*) RSTIMINC
      IF((NWS.EQ.1).AND.(NRS.EQ.1)) READ(15,*) RSTIMINC
      IF(ABS(NWS).EQ.2) THEN
        IF(NRS.EQ.0) READ(15,*) WTIMINC
        IF(NRS.EQ.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF(NWS.EQ.3) THEN
        READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC,&
     &             WREFTIM
        IF(NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,&
     &                          WLONINC,WTIMINC
        IF(NRS.EQ.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,&
     &                          WLONINC,WTIMINC,RSTIMINC
      ENDIF
      IF(ABS(NWS).EQ.4) THEN
        IF(NRS.EQ.0) READ(15,*) WTIMINC
        IF(NRS.EQ.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF(ABS(NWS).EQ.5) THEN
        IF(NRS.EQ.0) READ(15,*) WTIMINC
        IF(NRS.EQ.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF(NWS.EQ.6) THEN
        IF(NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,&
     &                          WLONINC,WTIMINC
        IF(NRS.EQ.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,&
     &                          WLONINC,WTIMINC,RSTIMINC
      ENDIF
      IF(NWS.EQ.10) THEN
        NWLAT=190
        NWLON=384
        IF(NRS.EQ.0) READ(15,*) WTIMINC
        IF(NRS.EQ.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF(NWS.EQ.11) THEN
        NWLAT=271
        NWLON=181
        WTIMINC=10800.
        IF(NRS.EQ.1) READ(15,*) RSTIMINC
       ENDIF
      READ(15,*) RNDAY
      READ(15,*) DRAMP
      READ(15,*) A00,B00,C00
      IF(NOLIFA.NE.2) THEN
        READ(15,*) H0
      ENDIF
      IF(NOLIFA.EQ.2) THEN
        READ(15,*) H0,NODEDRYMIN,NODEWETMIN,VELMIN
      ENDIF
      READ(15,*) SLAM0,SFEA0
!epm      SLAM0=SLAM0*DEG2RAD
!epm      SFEA0=SFEA0*DEG2RAD
      HBREAK=1.d0
      FTHETA=1.d0
      FGAMMA=1.d0
      IF(NOLIBF.EQ.0) READ(15,*) TAU
      CF=TAU
      IF(NOLIBF.EQ.1) READ(15,*) CF
      IF(NOLIBF.EQ.2) READ(15,*) CF,HBREAK,FTHETA,FGAMMA
      IF (IM.EQ.10) THEN
        READ(15,*) ESLM,ESLC
      ELSE
        READ(15,*) ESLM
      ENDIF
      READ(15,*) CORI
      READ(15,*) NTIF
      DO I=1,NTIF
!epm        READ(15,'(A5)')  TIPOTAG(I)
!epm        READ(15,*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
        READ(15,'(A5)')  
        READ(15,*)  
      END DO
      READ(15,*) NBFR
      DO I=1,NBFR
!epm        READ(15,'(A5)') BOUNTAG(I)
!epm        READ(15,*) AMIG(I),FF(I),FACE(I)
        READ(15,'(A5)') 
        READ(15,*) 
      END DO
      DO I=1,NBFR
        READ(15,'(A10)') ALPHA
        DO J=1,NETA
!epm          READ(15,*) EMO(I,J),EFA(I,J)
          READ(15,*) 
        END DO
      END DO
      READ(15,*) ANGINN
! epm assume NFLUXF=0
      NFLUXF = 0
      IF(NFLUXF.EQ.1) THEN
        READ(15,*) NFFR
      ENDIF
      READ(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      READ(15,*) NSTAE
      DO I=1,NSTAE
        IF(ICS.EQ.1) THEN
!epm          READ(15,*) XEL(I),YEL(I)
          READ(15,*) 
        ELSE
!epm          READ(15,*) SLEL(I),SFEL(I)
          READ(15,*) 
!epm          SLEL(I)=SLEL(I)*DEG2RAD
!epm          SFEL(I)=SFEL(I)*DEG2RAD
!epm          CALL CPP(XEL(I),YEL(I),SLEL(I),SFEL(I),SLAM0,SFEA0)
        ENDIF
      END DO
      READ(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      READ(15,*) NSTAV
      DO I=1,NSTAV
        IF(ICS.EQ.1) THEN
!epm          READ(15,*) XEV(I),YEV(I)
          READ(15,*) 
        ELSE
!epm          READ(15,*) SLEV(I),SFEV(I)
          READ(15,*) 
!epm          SLEV(I)=SLEV(I)*DEG2RAD
!epm          SFEV(I)=SFEV(I)*DEG2RAD
!epm          CALL CPP(XEV(I),YEV(I),SLEV(I),SFEV(I),SLAM0,SFEA0)
        ENDIF
      END DO
      IF(IM.EQ.10) THEN
        READ(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
        READ(15,*) NSTAC
        DO I=1,NSTAC
          IF(ICS.EQ.1) THEN
!epm            READ(15,*) XEC(I),YEC(I)
            READ(15,*) 
          ELSE
!epm            READ(15,*) SLEC(I),SFEC(I)
            READ(15,*) 
!epm            SLEC(I)=SLEC(I)*DEG2RAD
!epm            SFEC(I)=SFEC(I)*DEG2RAD
!epm            CALL CPP(XEC(I),YEC(I),SLEC(I),SFEC(I),SLAM0,SFEA0)
          ENDIF
        END DO
      ENDIF
      IF(NWS.NE.0) THEN
        READ(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
        READ(15,*) NSTAM
        DO I=1,NSTAM
          IF(ICS.EQ.1) THEN
!epm            READ(15,*) XEM(I),YEM(I)
            READ(15,*) 
          ELSE
!epm            READ(15,*) SLEM(I),SFEM(I)
            READ(15,*) 
!epm            SLEM(I)=SLEM(I)*DEG2RAD
!epm            SFEM(I)=SFEM(I)*DEG2RAD
!epm            CALL CPP(XEM(I),YEM(I),SLEM(I),SFEM(I),SLAM0,SFEA0)
          ENDIF
        END DO
      ENDIF
      READ(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      READ(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      IF(IM.EQ.10) THEN
        READ(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
      ENDIF
      IF(NWS.NE.0) THEN
        READ(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
      ENDIF
      READ(15,*) NFREQ 
      DO I=1,NFREQ  
!epm        READ(15,'(A10)') NAMEFR(I)
!epm        READ(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
        READ(15,'(A10)') 
        READ(15,*) 
      END DO
      READ(15,*) THAS,THAF,NHAINC,FMV
      READ(15,*) NHASE,NHASV,NHAGE,NHAGV
      READ(15,*) NHSTAR,NHSINC
      READ(15,*) ITITER,ISLDIA,CONVCR,ITMAX
      CLOSE(15)
      time_harmonic=THAF-THAS

      RETURN 
      END

