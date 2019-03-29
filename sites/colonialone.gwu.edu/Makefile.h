# FLASH makefile definitions for x86-64 Linux (intel compilers)
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

# Programming environment GNU or INTEL:
PE_ENV =GNU
#PE_ENV =INTEL


LIB_BASE   =/groups/balarasgrp/Software/$(PE_ENV)

MPI_PATH   =$(MPI_HOME)
HDF4_PATH  =
HDF5_PATH  =$(LIB_BASE)/hdf51.8.10

ZLIB_PATH  = 

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = 
MPE_PATH   =

#BLAS_PATH = $(LIB_BASE)/openblas
BLAS_PATH  = 

#LAPACK_PATH = $(LIB_BASE)/openblas
LAPACK_PATH  =

HYPRE_PATH   = $(LIB_BASE)/hypre2.9.0b
#HYPRE_PATH   = /groups/balarasgrp/akash/hypre
#HYPRE_PATH    = /home/akashdhruv/hypre-2.11.2/src/hypre

#SUPERLU_PATH = $(LIB_BASE)/SuperLU4.3
SUPERLU_PATH  = 

AMREX_PATH = /home/akashdhruv/amrex/amrex

ifeq      ($(NDIM), 2)
AMREX_PATH = /home/akashdhruv/amrex/amrex_2D
else ifeq ($(NDIM), 3)
AMREX_PATH = /home/akashdhruv/amrex/amrex_3D
else
AMREX_PATH = /home/akashdhruv/amrex/amrex_2D
endif

# Current directory:
export cur-dir := $(shell pwd)

# Set the location of top directory
export setup_dir = $(cur-dir)

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------

FCOMP   = ${MPI_PATH}/bin/mpif90
CCOMP   = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpicxx
LINK    = ${MPI_PATH}/bin/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for
#  flash_test, and is set for quick code generation, and (sometimes)
#  profiling.  The Makefile generated by setup will assign the generic token
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

ifeq ($(PE_ENV), INTEL)

# INTEL flags:
FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64 -xSSE4.2 -align array32byte
FFLAGS_DEBUG = -c -g -r8 -i4 -O0 -check bounds -check format \
-check output_conversion  -warn error -real_size 64 -check uninit \
-traceback -fp-stack-check -fpe0 -check pointers
FFLAGS_TEST  = ${FFLAGS_OPT} -fp-model precise
F90FLAGS = -DH5_USE_16_API -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE 
CFLAGS_DEBUG = -c -O0 -g -traceback -debug all -debug extended \
               -D_LARGEFILE64_SOURCE -ftrapuv -fp-stack-check
CFLAGS_TEST  = ${CFLAGS_OPT} -fp-model precise

else

# GNU Flags:
FFLAGS_OPT =  -c -O3 -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wuninitialized
FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fbacktrace -fdump-core -finit-real=nan \
-finit-integer=-999999 -fimplicit-none
FFLAGS_TEST =  -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none

CFLAGS_OPT   = -c -O3
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c

endif


#If we are using HDF5, we need to specify the path to the include files

FFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API

FLAGS_MPI   = -I$(MPI_PATH)/include

CFLAGS_BLAS = -I${BLAS_PATH}/include 
FFLAGS_BLAS = -I${BLAS_PATH}/include 

CFLAGS_LAPACK = -I${LAPACK_PATH}/include
FFLAGS_LAPACK = -I${LAPACK_PATH}/include

CFLAGS_SUPERLU = -I${SUPERLU_PATH}/include
FFLAGS_SUPERLU = -I${SUPERLU_PATH}/include

CFLAGS_HYPRE = -I${HYPRE_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

FFLAGS_AMREX = -I${AMREX_PATH}/include
FFLAGS_AMREX2D = ${FFLAGS_AMREX} -DN_DIM=2 -DNZB=1
#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

ifeq ($(PE_ENV), INTEL)

LFLAGS_OPT   = -xSSE4.2 -align array32byte -O3 -o
LFLAGS_DEBUG = -o
LFLAGS_TEST  = -O3 -o

else

LFLAGS_OPT   = -O3 -o
LFLAGS_DEBUG = -o
LFLAGS_TEST  = -O3 -o

endif


#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = 
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5 -lz

LIB_PAPI  =
ifeq ($(PE_ENV), INTEL)
LIB_MATH  = -limf -lm
else
LIB_MATH  = -lm
endif
LIB_MPI   = -lmpi_cxx 
LIB_NCMPI = -L${NCMPI_PATH}/lib -lpnetcdf
LIB_MPE   =

LIB_BLAS  = ${BLAS_PATH}/lib/libopenblas_sandybridgep-r0.2.3.a
LIB_LAPACK= ${BLAS_PATH}/lib/libopenblas_sandybridgep-r0.2.3.a
LIB_SUPERLU = -L${SUPERLU_PATH}/lib -lsuperlu_4.3 ${BLAS_PATH}/lib/libopenblas_sandybridgep-r0.2.3.a
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE  
LIB_STDCXX = -lstdc++
LIB_AMREX = -L${AMREX_PATH}/lib -lamrex
LIB_AMREX2D = ${LIB_AMREX}

#Specify TEC_PLOT=YES in order to link the tec plot library.
TEC_PLOT=NO
TEC_DIR=/groups/balarasgrp/Software/TecioLib
ifeq ($(TEC_PLOT), YES)
  CONFIG_LIB = -I${TEC_DIR} -L${TEC_DIR}  -ltecio 
endif

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ =

#----------------------------------------------------------------------------
# Additional commands
#----------------------------------------------------------------------------

MV = mv -f
AR = xiar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo

ifeq ($(FLASHBINARY),true)


FFLAGS_WO_WARNALL = $(patsubst -pedantic,,$(FFLAGS))

#Files mix and match assumed shape arrays, assumed size arrays
#and scalars in function calls.  This is fine but it is viewed as
#a problem when using strict type checking compiler options.
fftpack.o : %.o : %.f90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<
gr_pfftDcftForward.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<
gr_pfftDcftInverse.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<

endif


#Configure lines:
#
# mpich2-1.2.1p1:
#./configure --prefix=/opt/mpich2/1.2.1p1/gcc-4.4.3
#--enable-error-checking=all --enable-error-messages=all 
#--with-pm=gforker:mpd --enable-g=dbg,meminit --enable-fast=defopt 
#--enable-f77 --enable-f90 --enable-cxx --enable-romio --enable-sharedlibs=gcc 
#--with-mpe CC=gcc F77=gfortran F90=gfortran CXX=g++ 2>&1 | tee ../mpich2_1.2.1p1_gcc-4.4.3_build.out
#
# hdf5-1.8.4-patch1:
#./configure --prefix=/opt/hdf5/1.8.4-patch1/gcc-4.4.3 
#CC=mpicc FC=mpif90 CXX=mpicxx --enable-production --enable-debug=all --enable-shared
#--enable-parallel --enable-using-memchecker 2>&1 | tee ../hdf5_1.8.4-patch1_gcc-4.4.3_build.out
#
# parallel-netcdf-1.1.1:
#./configure --prefix=/opt/parallel-netcdf/1.1.1/gcc-4.4.3 
#--enable-fortran --with-mpi=/opt/mpich2/1.2.1p1/gcc-4.4.3 
#CFLAGS="${CFLAGS} -g" 2>&1 | tee ../parallel-netcdf-1.1.1_gcc-4.4.3_build.out
#
# valgrind-3.5.0: (Need to patch because we have a new version of glibc)
#patch -Np0 -i ../valgrind_glibc211.diff || return 1
#autoreconf
#./configure --prefix=/opt/valgrind/3.5.0/gcc-4.4.3
#--without-mpicc 2>&1 | tee ../valgrind-3.5.0_gcc-4.4.3_build.out