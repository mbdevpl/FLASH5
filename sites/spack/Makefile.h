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
LINK    = mpif90 -std=c++11
# -fopenmp

# pre-processor flag
PP      = -D


IS_FCOMP_GNU := $(shell ${FCOMP} --version | grep "GNU\\|GCC" | wc -l)
IS_CCOMP_GNU := $(shell ${CCOMP} --version | grep "GNU\\|GCC" | wc -l)
# @echo Is Fortran compiler GNU?
# @echo $(value IS_FCOMP_GNU)
$(info Is Fortran compiler GNU? $(value IS_FCOMP_GNU))
$(info Is C compiler GNU? $(value IS_CCOMP_GNU))

IS_FCOMP_PGI := $(shell ${FCOMP} --version | grep "PGI" | wc -l)
IS_CCOMP_PGI := $(shell ${CCOMP} --version | grep "PGI" | wc -l)
$(info Is Fortran compiler PGI? $(value IS_FCOMP_PGI))
$(info Is C compiler PGI? $(value IS_CCOMP_PGI))


ifneq (${IS_FCOMP_GNU}, 0)

ifneq (${IS_FCOMP_PGI}, 0)
$(error Fortran compiler detected as both GNU and PGI!)
endif

FFLAGS_OPT = -g -c -O2 -fdefault-real-8 -fdefault-double-8 \
$(shell echo -I${CPATH} | sed 's|:| -I|g') \
-Wuninitialized
# -fopenmp \

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

else
ifneq (${IS_FCOMP_PGI}, 0)

LINK = mpif90 \
-fast -Mvect=simd -Mcache_align -Mflushz -Mpre \
-acc -ta=tesla:cc60 -ta=tesla:nordc -Minfo=all
# -mp

FFLAGS_OPT = -g -c -O 4 -r8 -i4 \
$(shell echo -I${CPATH} | sed 's|:| -I|g') \
-fast -Mvect=simd -Mcache_align -Mflushz -Mpre \
-acc -ta=tesla:cc60 -ta=tesla:nordc -Minfo=all
# -mp

FFLAGS_DEBUG = -g -c -O0 -r8 -i4 \
$(shell echo -I${CPATH} | sed 's|:| -I|g')

FFLAGS_TEST = -g -c -r8 -i4 \
$(shell echo -I${CPATH} | sed 's|:| -I|g')

else
$(error Fortran compiler detected as neither GNU nor PGI!)
endif
endif

FFLAGS_HYPRE =
FFLAGS_AMREX =
FFLAGS_AMREX2D = -DN_DIM=2 -DNZB=1

F90FLAGS =


ifneq (${IS_CCOMP_GNU}, 0)

ifneq (${IS_CCOMP_PGI}, 0)
$(error C compiler detected as both GNU and PGI!)
endif

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

else
ifneq (${IS_CCOMP_PGI}, 0)

CFLAGS_OPT = -g -c -O2

CFLAGS_DEBUG = -g -c -O0

CFLAGS_TEST = -g -c

else
$(error C compiler detected as neither GNU nor PGI!)
endif
endif



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

LIB_MPI   = -lmpi
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
