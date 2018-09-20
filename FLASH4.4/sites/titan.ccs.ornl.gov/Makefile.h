# FLASH makefile definitions for Cray XT4 System, Jaguar at ORNL
#
# The XT4 makes use of modules.  
# Load the following modules before compiling (names as of 5th Nov 08).
#
# HDF5: module load hdf5/1.6.7_par
# PNETCDF: module load p-netcdf/1.0.2 
# TAU: module load tau/2.17.2 
#     use -tau=${TAUROOT}/lib/Makefile.tau-pgi-callpath-mpi-pdt at setup time.

#----------------------------------------------------------------------------
# Set the HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

LC_PE_ENV = $(shell echo ${PE_ENV} | tr A-Z a-z)

MPI_PATH   =
PAPI_PATH  = 
PAPI_FLAGS = 
NCMPI_PATH = 
AMREX_PATH = ${AMREX_DIR}/titan_${NDIM}d_${LC_PE_ENV}

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifeq ($(findstring $(PE_ENV),INTEL GNU PGI CRAY),)
$(error Your environment "$(PE_ENV)" is invalid.  It must be "INTEL", "GNU", "PGI" or "CRAY")
else
$(warning You are using the "$(PE_ENV)" environment)
endif

ifdef CRAY_HDF5_DIR
$(warning You are using the HDF5 installation at "$(CRAY_HDF5_DIR)")
else
$(warning Generally FLASH needs HDF5.  You can load it using "module load cray-hdf5-parallel")
endif


FCOMP   =  ftn
CCOMP   =  cc
CPPCOMP =  CC
CUCOMP  =  nvcc
LINK    =  ftn

#-----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized code
#  code ("-opt"), one for debugging ("-debug"), and one for testing ("-test").
#  Passing these flags to the setup script will cause the value associated with
#  the corresponding keys (i.e. those ending in "_OPT", "_DEBUG", or "_TEST") to
#  be incorporated into the final Makefile. For example, passing "-opt" to the
#  setup script will cause the flags following "FFLAGS_OPT" to be assigned to
#  "FFLAGS" in the final Makefile. If none of these flags are passed, the default
#  behavior will match that of the "-opt" flag.
#  In general, "-opt" is meant to optimize compilation and linking. "-debug"
#  should enable runtime bounds checking, debugger symbols, and other compiler-
#  specific debugging options. "-test" is useful for testing different
#  combinations of compiler flags particular to your individual system.
#----------------------------------------------------------------------------

ifeq ($(PE_ENV),CRAY)

    OPENMP       =

    FFLAGS_OPT   = -c -s real64 -s integer32 -O2 -eZ
    FFLAGS_DEBUG = -c -s real64 -s integer32 -O0 -eD -eZ
    FFLAGS_TEST  = -c -s real64 -s integer 32 -em

    CFLAGS_OPT   = -c -O2
    CFLAGS_DEBUG = -c -O0 -eD
    CFLAGS_TEST  = -c 

    LFLAGS_OPT   = -O2 -o
    LFLAGS_DEBUG = -eD -o
    LFLAGS_TEST  = -em -o

    FFLAGS_OACC  =
    CFLAGS_OACC  =
    LIB_OACC     =

else ifeq ($(PE_ENV),INTEL)

    OPENMP       = -qopenmp

    FFLAGS_OPT   = -c -real-size 64 -integer-size 32 -O2 -fp-model precise -fpp
#This is an alternative set of debug flags that doesn't do any checking -- a lot of the
#checks fail and aren't useful if you just want symbols for debugging
#FFLAGS_DEBUG = -c -real-size 64 -integer-size 32 -O0 -g -fpp
    FFLAGS_DEBUG = -c -real-size 64 -integer-size 32 -O0 -g -fpp -traceback -check uninit -check bounds -ftrapuv
    FFLAGS_TEST  = -c -real-size 64 -integer-size 32 -list

    CFLAGS_OPT   = -c -O2
    CFLAGS_DEBUG = -c -g -O0
    CFLAGS_TEST  = -c 

    LFLAGS_OPT   = -O2 -fp-model precise -Wl,--whole-archive,-ldmapp,--no-whole-archive -o
    LFLAGS_DEBUG = -g -o
    LFLAGS_TEST  = -list -o

    FFLAGS_OACC  =
    CFLAGS_OACC  =
    LIB_OACC     =

else ifeq ($(PE_ENV),PGI)

    OPENMP       = -mp

    FFLAGS_OPT   = -c -r8 -i4 -fastsse -tp bulldozer -Mpreprocess# -DDEBUG_DRIVER
    FFLAGS_DEBUG = -c -r8 -i4 -O0 -g -Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mdclchk -Mbounds -Mpreprocess
    FFLAGS_TEST  = -c -r8 -i4 -Mprof=lines

    CFLAGS_OPT   = -c -fastsse -tp bulldozer
    CFLAGS_DEBUG = -c -O0 -g
    CFLAGS_TEST  = -c 

    LFLAGS_OPT   = -fastsse -tp bulldozer -Wl,--whole-archive,-ldmapp,--no-whole-archive -o
    LFLAGS_DEBUG = -g -o
    LFLAGS_TEST  = -Mprof=lines -o

    FFLAGS_OACC  = -acc -ta=tesla:cc35
    CFLAGS_OACC  = -acc -ta=tesla:cc35
    LIB_OACC     = -acc -ta=tesla:cc35 -acclibs

else ifeq ($(PE_ENV),GNU)

    OPENMP       = -fopenmp

    FFLAGS_OPT   = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -O2 -cpp
    FFLAGS_DEBUG = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -Og -g -fbacktrace -fcheck=bounds,do,mem,pointer -cpp
    FFLAGS_TEST  = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

    CFLAGS_OPT   = -c -O2
    CFLAGS_DEBUG = -c -Og -g -fbounds-check
    CFLAGS_TEST  = -c 

    LFLAGS_OPT   = -O2 -o
    LFLAGS_DEBUG = -g -o
    LFLAGS_TEST  = -o

    FFLAGS_OACC  =
    CFLAGS_OACC  =
    LIB_OACC     =

endif

CU_FLAGS = -c -g -O2 -m64 -gencode arch=compute_35,code=sm_35

#No path required because we are using compiler wrapper scripts.
FFLAGS_MPI   =
CFLAGS_MPI   =

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = ${HDF5_CLIB} -DH5_USE_16_API
CFLAGS_NCMPI = ${PNETCDF_LIB}

FFLAGS_PAPI  = 

FFLAGS_AMREX = -I${AMREX_PATH}/include

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

LIB_HDF5    = ${HDF5_CLIB}
LIB_LAPACK  =
LIB_CUDA    = -lcublas -lcudart -lcuda

LIB_AMREX   = -L${AMREX_PATH}/lib -lamrex -lstdc++
LIB_STDCXX = -lstdc++
              
LIB_MPI     = 
LIB_PAPI    = 
LIB_PNG     = 

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   = ${PNETCDF_LIB}
LIB_MPE     =
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

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo
AWK   = awk
CAT   = cat
