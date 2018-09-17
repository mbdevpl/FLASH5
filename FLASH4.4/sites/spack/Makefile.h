# FLASH makefile definitions avoiding any hardcoded paths.
# This is useful when using Spack, and might also be useful in custom environments.

MPI_PATH   =
HDF5_PATH  =
HYPRE_PATH =
AMREX_PATH =
ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

LIB_NCMPI = /usr/local
MPE_PATH   =

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC
LINK    = mpif90 -std=c++11 -fopenmp

# pre-processor flag
PP      = -D



FFLAGS_OPT = -g -c -O2 -fdefault-real-8 -fdefault-double-8 \
$(shell echo -I${CPATH} | sed 's|:| -I|g') \
-Wuninitialized

FFLAGS_DEBUG = -g -ggdb -c -O0 -fdefault-real-8 -fdefault-double-8 \
$(shell echo -I${CPATH} | sed 's|:| -I|g') \
-pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fimplicit-none -fstack-protector-all \
-fbacktrace -fbounds-check

FFLAGS_TEST = -g -ggdb -c -fdefault-real-8 -fdefault-double-8 \
$(shell echo -I${CPATH} | sed 's|:| -I|g') \
-ffree-line-length-none

FFLAGS_HYPRE =
FFLAGS_AMREX =
FFLAGS_AMREX2D = -DN_DIM=2 -DNZB=1


F90FLAGS =


#The macro _FORTIFY_SOURCE adds some lightweight checks for buffer
#overflows at both compile time and run time (only active at -O1 or higher)
#http://gcc.gnu.org/ml/gcc-patches/2004-09/msg02055.html
CFLAGS_OPT = -g -c -O2 \
-Wuninitialized -D_FORTIFY_SOURCE=2

CFLAGS_DEBUG = -g -ggdb -c -O0 \
-Wno-div-by-zero -Wundef \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-protector-all

CFLAGS_TEST = -g -c


# Platform symbol
CDEFINES += -DDarwin

CFLAGS_HDF5 = -DH5_USE_16_API
CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -g \
$(shell echo -L${LD_LIBRARY_PATH} | sed 's|:| -L|g') \
-o
LFLAGS_DEBUG = -g -O0 \
$(shell echo -L${LD_LIBRARY_PATH} | sed 's|:| -L|g') \
-o
LFLAGS_TEST  = -g \
$(shell echo -L${LD_LIBRARY_PATH} | sed 's|:| -L|g') \
-o


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

LIB_HDF5  = -lhdf5 -lhdf5_fortran

LIB_BLAS   = -lopenblas
LIB_SUPERLU= -lsuperlu
LIB_LAPACK = $(LIB_BLAS) $(LIB_SUPERLU)

LIB_PAPI  =
LIB_MATH  =

LIB_MPI   =
LIB_MPE   =

LIB_HYPRE = -lHYPRE

LIB_AMREX = -lamrex
LIB_AMREX2D = ${LIB_AMREX}
LIB_STDCXX =
LIB_STDCXX = -lstdc++

# Uncomment the following line to use electic fence memory debugger.
# export EF_ALLOW_MALLOC_0=1
#CONFIG_LIB = -L/usr/lib64 -lefence

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
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo


#----------------------------------------------------------------------------
# Fake existence of iso_c_bindings module to prevent unnecessary recompilations.
#----------------------------------------------------------------------------
#ifeq ($(FLASHBINARY),true)
#iso_c_binding.mod :
#	touch $@
#endif
