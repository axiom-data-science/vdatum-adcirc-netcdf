#***********************************************************************
# Makefile   
#************ ***********************************************************
# This makefile may be used to create an executable test mainline
# for the Hydro_NetCDF library
#
# make maintestfem
# make maintestgrid
# make clean
#
# ?= assignments allow for command line overrides of makefile variables
# tfg  didn't get the ?= to work on my RH9,0
#
# 2004-08-27:
# built on Fedora FC2 & Debian testing with intel ifort 8.0
# FC='ifort -static-libcxa' F77=$(FC) FFLAGS='-I/usr/local/include' make 
#   all works fine -- Dave Forrest <drf5n@vims.edu>
 
FFLAGS = -O

EXEC = adcirc2netcdf
AD2CDF = adcirc2netcdf

CDFLIBDIR = $(CONDA_PREFIX)/lib
FC = $(CONDA_PREFIX)/bin/gfortran
CDFINCL = -I$(CONDA_PREFIX)/include
CDFLIBS = -L$(CDFLIBDIR) -lnetcdf -lhdf5 -lhdf5_hl -lnetcdff
 
SRCS = tides_netcdf.f adcirc2netcdf.f90

OBJS = ${SRCS:.f,.f90=.o}

$(AD2CDF):	$(OBJS)
	$(FC) -o $(AD2CDF) $(FFLAGS) $(OBJS) $(CDFINCL) $(CDFLIBS) 

convert: $(AD2CDF)
	LD_LIBRARY_PATH=$(CDFLIBDIR) ./adcirc2netcdf
