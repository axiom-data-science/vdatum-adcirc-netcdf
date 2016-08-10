ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tidescdf_fem.f  Writes standardized NetCDF files for tides  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write_tidescdf_fem(netcdf_file,ncid,globalstr,
     & ne,nn,nvrt,ntides,ele,lon,lat,depth,tidenames,tidefreqs,
     & ampl,phase,u_amp,u_phase,v_amp,v_phase,
     & tideanalysis,timeanalysis)

      implicit none

      include 'netcdf.inc'
c input variables
      integer ne,nn,nface,ntides,nvrt
      integer ele(3,ne),j
      real lon(nn),lat(nn),depth(nn),
     & ampl(ntides,nn),phase(ntides,nn),timeanalysis,tidefreqs(ntides)
      real u_amp(ntides,nn,nvrt),v_amp(ntides,nn,nvrt)
      real u_phase(ntides,nn,nvrt),v_phase(ntides,nn,nvrt)
      character globalstr(9)*40
      character tideanalysis*20
      character tidenames(ntides)*10
      character netcdf_file*80

C Netcdf internal variables
      integer iret,ncid,intval(4),CORNER(4),COUNT(4)
      integer ele_id, nele_dim, char_dim
      integer node_dim, nface_dim, ntides_dim, nvrt_dim
      integer lon_id, lat_id, depth_id, nvrt_id
      integer tidefreqs_id, ampl_id, phase_id
      integer uamp_id, uphase_id, vamp_id, vphase_id
      integer tidenames_id
      integer values(8)

      logical lu_amp,lu_phase,lv_amp,lv_phase

C Set optional variable flags
      lu_amp = .TRUE.
      lu_phase = .TRUE.
      lv_amp = .TRUE.
      lv_phase = .TRUE.
      if (u_amp(1,1,1).le.0.0) lu_amp=.FALSE.
      if (u_phase(1,1,1).le.0.0) lu_phase=.FALSE.
      if (v_amp(1,1,1).le.0.0) lv_amp=.FALSE.
      if (v_phase(1,1,1).le.0.0) lv_phase=.FALSE.
      
c      save  ele_id, nele_dim
c      save  lon_id, lat_id, depth_id
c      save  ampl_id, phase_id
c      save  tidenames_id, tidefreqs_id
      
C Check for Triangular mesh or Quadralateral or Mixed
      if (globalstr(1) .eq. "Triangular") then
          nface=3
      else
          nface=4
      endif

C Initialize
c open file
      iret = nf_create(netcdf_file, NF_CLOBBER, ncid)
      call check_err(iret)
c define dimensions
      iret = nf_def_dim(ncid, 'node', nn, node_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'nele', ne, nele_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'nface', nface, nface_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'ntides', ntides, ntides_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'charlength', 10, char_dim)
      call check_err(iret)
      if (lu_amp .or. lu_phase .or. 
     &    lv_amp .or. lv_phase) then
        if (globalstr(2) .eq. "Z" ) then
          iret = nf_def_dim(ncid, 'zgrid', nvrt, nvrt_dim)
          call check_err(iret)
        elseif (globalstr(2) .eq. "2D") then
          iret = nf_def_dim(ncid, 'depth-averaged', nvrt, nvrt_dim)
          call check_err(iret)
        else
          iret = nf_def_dim(ncid, 'sigma', nvrt, nvrt_dim)
          call check_err(iret)
        endif
      endif


c define variables
C tidenames
      intval(2) = ntides_dim
      intval(1) = char_dim
      iret = nf_def_var(ncid, 'tidenames',NF_CHAR,2,intval,
     & tidenames_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, tidenames_id, 'long_name'
     & , 16, 'Tide Constituent')
      iret = nf_put_att_int(ncid, tidenames_id, 
     & 'missing_value', NF_INT, 1,-1)
      iret = nf_put_att_text(ncid, tidenames_id, 'standard_name', 16, 
     & 'tide_constituent')
      call check_err(iret)
C tidefreqs
      iret = nf_def_var(ncid, 'tidefreqs', NF_REAL,
     & 1, ntides_dim, tidefreqs_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, tidefreqs_id, 'long_name'
     & , 14, 'Tide Frequency')
      iret = nf_put_att_text(ncid, tidefreqs_id, 'units', 14, 
     & 'radians/second')
      call check_err(iret)
      iret = nf_put_att_real(ncid, tidefreqs_id, 
     & 'missing_value', NF_REAL, 1, -99999.0)
      iret = nf_put_att_text(ncid, tidefreqs_id, 'standard_name', 14, 
     & 'tide_frequency')
      call check_err(iret)
C ele
      intval(1) = nface_dim
      intval(2) = nele_dim
      iret = nf_def_var(ncid, 'ele', NF_INT, 2, intval, ele_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, ele_id, 'long_name'
     & , 33,'Horizontal Element Incidence List')
      iret = nf_put_att_int(ncid, ele_id, 
     & 'missing_value', NF_INT, 1,-1)
C lon
      intval(1) = node_dim
      iret = nf_def_var(ncid, 'lon', NF_REAL, 1, intval, lon_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9, 'Longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'standard_name', 9, 
     & 'longitude')
C lat
      iret = nf_def_var(ncid, 'lat', NF_REAL, 1, intval, lat_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8, 'Latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lat_id, 'standard_name', 8, 
     & 'latitude')
C depth
      intval(1) = node_dim
      iret = nf_def_var(ncid, 'depth', NF_REAL, 1, intval, depth_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid,depth_id,'long_name',11,'Bathymetry')
      call check_err(iret)
      iret = nf_put_att_text(ncid, depth_id,'units',6,'meters')
      iret = nf_put_att_text(ncid, depth_id,'positive',4,'down')
      call check_err(iret)
       iret = nf_put_att_text(ncid, depth_id, 'standard_name', 5, 
     & 'depth')

      if (lu_amp .or. lu_phase .or. 
     &    lv_amp .or. lv_phase) then
C alternative attributes for sigma or zlevel
      if (globalstr(2) .eq. "Z" ) then
C Zgrid
        intval(1) = nvrt_dim
        iret = nf_def_var(ncid, 'zgrid', NF_REAL, 1, intval, nvrt_id)
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'long_name', 25,
     &   'Z Coordinate Fixed Levels')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'units', 1, 'm')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'positive', 2,'up')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'standard_name', 22,
     &   'ocean_z_coordinate')
      elseif (globalstr(2) .eq. "2D" ) then
C Depth-Averaged
        intval(1) = nvrt_dim
        iret = nf_def_var(ncid, 'depth-averaged', NF_REAL, 1, 
     &   intval, nvrt_id)
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'long_name', 37,
     &   'Depth-Averaged Output in the Vertical')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'standard_name', 37,
     &   'depth_averaged_output_in_the_vertical')
      else
C sigma
        intval(1) = nvrt_dim
        iret = nf_def_var(ncid, 'sigma', NF_REAL, 1,intval, nvrt_id)
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'long_name', 44,
     &   'Sigma Stretched Vertical Coordinate at Nodes')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'units', 11, 
     &   'sigma_level')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'positive', 4, 'down')
        call check_err(iret)
        iret = nf_put_att_text(ncid, nvrt_id, 'standard_name', 22,
     &   'ocean_sigma_coordinate')
        iret = nf_put_att_text(ncid, nvrt_id, 'formula_terms', 35,
     &   'sigma: sigma eta: zeta depth: depth')
       endif
       endif

C tide amplitude
      intval(1) = ntides_dim
      intval(2) = node_dim
      iret = nf_def_var(ncid, 'ampl', NF_REAL, 2, intval, ampl_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, ampl_id, 'long_name', 14, 
     & 'Tide Amplitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ampl_id, 'units', 6, 'meters')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ampl_id, 'tide_analysis', 20, 
     & tideanalysis)
      call check_err(iret)
      iret = nf_put_att_real(ncid, ampl_id, 'time_analysis', 1, 
     & timeanalysis)
      iret = nf_put_att_real(ncid, ampl_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, ampl_id, 
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, ampl_id, 'standard_name', 14, 
     & 'tide_amplitude')

C tide phase
      intval(1) = ntides_dim
      intval(2) = node_dim
      iret = nf_def_var(ncid, 'phase', NF_REAL, 2,intval, phase_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, phase_id, 'long_name', 10, 
     & 'Tide Phase')
      call check_err(iret)
      iret = nf_put_att_text(ncid, phase_id, 'units', 11, 'degrees UTC')
      call check_err(iret)
      iret = nf_put_att_text(ncid, phase_id, 'tide_analysis', 20, 
     & tideanalysis)
      call check_err(iret)
      iret = nf_put_att_real(ncid, phase_id, 'time_analysis', 1, 
     & timeanalysis)
      iret = nf_put_att_real(ncid, phase_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, phase_id, 
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, phase_id, 'standard_name', 10, 
     & 'tide_phase')

C u amplitude
      if (lu_amp) then
      intval(3) = nvrt_dim
      intval(2) = node_dim
      intval(1) = ntides_dim
      iret = nf_def_var(ncid, 'u_amp', NF_REAL, 3, intval, uamp_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, uamp_id, 'long_name'
     & , 33,'Eastward Water Velocity Amplitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, uamp_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, uamp_id,
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, uamp_id,
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, uamp_id, 'standard_name', 37,
     & 'eastward_sea_water_velocity_amplitude')
      endif

C v amplitude
      if (lv_amp) then
C     intval(1) = ntides_dim
C     intval(2) = node_dim
C     intval(3) = nvrt_dim
      intval(2) = nvrt_dim
      intval(2) = node_dim
      intval(1) = ntides_dim
      iret = nf_def_var(ncid, 'v_amp', NF_REAL, 3, intval, vamp_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, vamp_id, 'long_name'
     & , 34,'Northward Water Velocity Amplitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, vamp_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, vamp_id,
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, vamp_id,
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, vamp_id, 'standard_name', 38,
     & 'northward_sea_water_velocity_amplitude')
      endif

C u phase
      if (lu_phase) then
C     intval(1) = ntides_dim
C     intval(2) = node_dim
C     intval(3) = nvrt_dim
      intval(3) = nvrt_dim
      intval(2) = node_dim
      intval(1) = ntides_dim
      iret = nf_def_var(ncid, 'u_phase', NF_REAL, 3, intval, uphase_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, uphase_id, 'long_name'
     & , 29,'Eastward Water Velocity Phase')
      call check_err(iret)
      iret = nf_put_att_text(ncid, uphase_id, 'units', 3, 'degrees UTC')
      call check_err(iret)
      iret = nf_put_att_real(ncid, uphase_id,
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, uphase_id,
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, uphase_id, 'standard_name', 33,
     & 'eastward_sea_water_velocity_phase')
      endif

C v phase
      if (lv_phase) then
C     intval(1) = ntides_dim
C     intval(2) = node_dim
C     intval(3) = nvrt_dim
      intval(3) = nvrt_dim
      intval(2) = node_dim
      intval(1) = ntides_dim
      iret = nf_def_var(ncid, 'v_phase', NF_REAL, 3, intval, vphase_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, vphase_id, 'long_name'
     & , 30,'Northward Water Velocity Phase')
      call check_err(iret)
      iret = nf_put_att_text(ncid, vphase_id, 'units', 3, 'degrees UTC')
      call check_err(iret)
      iret = nf_put_att_real(ncid, vphase_id,
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_real(ncid, vphase_id,
     & '_FillValue', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, vphase_id, 'standard_name', 34,
     & 'northward_sea_water_velocity_phase')
      endif

C Global Attributes
      iret = nf_put_att_text(ncid, NF_GLOBAL,'file_type'
     & ,3,'FEM')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'Conventions'
     & ,6,'CF-1.0')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'grid_type'
     & ,len_trim(globalstr(1)),globalstr(1))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'z_type'
     & ,len_trim(globalstr(2)),globalstr(2))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'model'
     & ,len_trim(globalstr(3)),globalstr(3))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'title'
     & ,len_trim(globalstr(4)),globalstr(4))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'comment'
     & ,len_trim(globalstr(5)),globalstr(5))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'source'
     & ,len_trim(globalstr(6)),globalstr(6))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'institution'
     & ,len_trim(globalstr(7)),globalstr(7))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'history'
     & ,len_trim(globalstr(8)),globalstr(8))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'references'
     & ,len_trim(globalstr(9)),globalstr(9))
      call check_err(iret)
      
c  END DEFINITIONS 
      iret = nf_enddef(ncid)
      call check_err(iret)
      write(*,*) 'End definitions'

c  write lon,lat,depth,ele
      CORNER(1)=1
      COUNT(1)=nn
      iret=nf_put_vara_real(ncid,lon_id,CORNER,COUNT,lon)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,lat_id,CORNER,COUNT,lat)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,depth_id,CORNER,COUNT,depth)
      call check_err(iret)

      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=nface
      COUNT(2)=ne
      iret=nf_put_vara_int(ncid,ele_id,CORNER,COUNT,ele)
      write(*,*) ' done ele'
      call check_err(iret)

      write(*,*) 'Initialize done',COUNT(1),COUNT(2)

c  write tidenames, tidefreqs
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=10
      COUNT(2)=ntides
      iret=nf_put_vara_text(ncid,tidenames_id,CORNER,COUNT,tidenames)
      call check_err(iret)
      CORNER(1) = 1
      COUNT(1)=ntides
      iret=nf_put_vara_real(ncid,tidefreqs_id,CORNER,COUNT,tidefreqs)
      call check_err(iret)

c  write ampl, phase
c  2-D fields nodal
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(2)=nn
      COUNT(1)=ntides
      iret=nf_put_vara_real(ncid,ampl_id,CORNER,COUNT,ampl)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,phase_id,CORNER,COUNT,phase)
      call check_err(iret)
      write(*,*) 'starting the velocity write',nvrt,nn,ntides

c  write ampl, phase for velocities
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = 1
      COUNT(3)=nvrt
      COUNT(2)=nn
      COUNT(1)=ntides
      if (lu_amp) then
      iret=nf_put_vara_real(ncid,uamp_id,CORNER,COUNT,u_amp)
      call check_err(iret)
      endif
      if (lu_phase) then
      iret=nf_put_vara_real(ncid,uphase_id,CORNER,COUNT,u_phase)
      call check_err(iret)
      endif
      if (lv_amp) then
      iret=nf_put_vara_real(ncid,vamp_id,CORNER,COUNT,v_amp)
      call check_err(iret)
      endif
      if (lv_phase) then
      iret=nf_put_vara_real(ncid,vphase_id,CORNER,COUNT,v_phase)
      call check_err(iret)
      endif

c  close file
      iret = nf_close(ncid)
      call check_err(iret)
      write(*,*) 'Close ',netcdf_file 
      
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
