# /homes/joneal/.soft.cache.sh
# Created on Fri Mar  9 14:42:25 2018.
#
# DO NOT MODIFY THIS FILE DIRECTLY!!
# If you do, your modifications will probably be overwritten.
#
# This file is created automatically by the soft system.  It is used
# to build your environment.  It was generated from the contents of 
# this file:
#
#  /homes/joneal/.soft
#
# If the contents of this file aren't correct, then modify that file, and 
# type 'resoft'.  The contents of this file are listed below, as are
# any errors found when parsing it.
#
# More information about the SoftEnv system can be found by typing 
# 'man softenv-intro'.
#
# =========================================================================
#
# Contents of /homes/joneal/.soft:
#
# #
# # This is your SoftEnv configuration run control file.
# #
# #   It is used to tell SoftEnv how to customize your environment by
# #   setting up variables such as PATH and MANPATH.  To learn more
# #   about this file, do a "man softenv".
# #
# 
# #+git-2.10.1
# +python-2.7.9
# #PYTHONPATH=$HOME/local/lib/python2.7/site-packages
# #+cmake-3.5.1
# +szip-2.1-gcc-6.2.0
# 
# #+intel-17-initial
# #+mkl-17-initial
# #+mpich-intel
# 
# +gcc-6.2.0
# +mpich-3.2-gcc-6.2.0
# +hdf5-1.8.20-gcc-6.2.0-mpich-3.2-parallel-fortran
# 
# #-------------------------------------------------------------------------------
# @default
# #-------------------------------------------------------------------------------
# #+ddt
# #+tau-2.25.1
# #+totalview
# #+texlive-2015
# #+matlab-2016b
# 
#
# ===========================================================================
#
# Errors while reading this file:
#
# None.#
# ===========================================================================

# Initialize CPLUS_INCLUDE_PATH
CPLUS_INCLUDE_PATH=''
export CPLUS_INCLUDE_PATH

# Initialize C_INCLUDE_PATH
C_INCLUDE_PATH=''
export C_INCLUDE_PATH

# Initialize IDL_PATH
IDL_PATH=''
export IDL_PATH

# Initialize INTEL_LICENSE_FILE
INTEL_LICENSE_FILE=''
export INTEL_LICENSE_FILE

# Initialize LD_LIBRARYN32_PATH
LD_LIBRARYN32_PATH=''
export LD_LIBRARYN32_PATH

# Initialize LD_LIBRARYN64_PATH
LD_LIBRARYN64_PATH=''
export LD_LIBRARYN64_PATH

# Initialize LD_LIBRARY_PATH
LD_LIBRARY_PATH=''
export LD_LIBRARY_PATH

# Initialize LM_LICENSE_FILE
LM_LICENSE_FILE=''
export LM_LICENSE_FILE

# Initialize MANPATH
MANPATH=''
export MANPATH

# Initialize MATLABPATH
MATLABPATH=''
export MATLABPATH

# Initialize MKL_HOME
MKL_HOME=''
export MKL_HOME

# Initialize NAG_KUSARI_FILE
NAG_KUSARI_FILE=''
export NAG_KUSARI_FILE

# Initialize PATH
PATH=''
export PATH

# Initialize PYTHONPATH
PYTHONPATH=''
export PYTHONPATH

# Initialize SIDL_DLL_PATH
SIDL_DLL_PATH=''
export SIDL_DLL_PATH

# Initialize TEXINPUTS
TEXINPUTS=''
export TEXINPUTS

# translating: initvar PATH
PATH=/INIT:/bin:/sbin:/usr/local/bin

# translating: initvar WHATAMI
WHATAMI=`/bin/whatami`

# translating: initvar ARCH
ARCH=`/bin/whatami`

# translating: switch ${ARCH}
case ${ARCH} in
  
  irix-6)
    
    # translating: initvar IRIX_ARCH
    IRIX_ARCH=`uname`
  
  ;;
  #  End of irix-6 case choice
  
  solaris-2)
    
    # translating: initvar SOLARIS_ARCH
    SOLARIS_ARCH=`uname -m`
    
    # translating: initvar SOLARIS_VERSION
    SOLARIS_VERSION=`uname -r`
  
  ;;
  #  End of solaris-2 case choice
  
  solaris-9)
    
    # translating: initvar SOLARIS_ARCH
    SOLARIS_ARCH=`uname -m`
    
    # translating: initvar SOLARIS_VERSION
    SOLARIS_VERSION=`uname -r`
  
  ;;
  #  End of solaris-9 case choice

esac
# End of case $${ARCH}

# translating: switch ${ARCH}
case ${ARCH} in
  
  aix-4)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of aix-4 case choice
  
  aix-5)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of aix-5 case choice
  
  irix-5)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of irix-5 case choice
  
  irix-6)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of irix-6 case choice
  
  linux)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of linux case choice
  
  linux-2)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-2 case choice
  
  linux-Ubuntu_10.04-ia32)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-Ubuntu_10.04-ia32 case choice
  
  linux-Ubuntu_10.04-x86_64)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-Ubuntu_10.04-x86_64 case choice
  
  linux-Ubuntu_12.04-ia32)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_12.04-ia32 case choice
  
  linux-Ubuntu_12.04-x86_64)
    
    # translating: append LD_LIBRARY_PATH
    if [ "x${LD_LIBRARY_PATH}" = "x" ]; then
      LD_LIBRARY_PATH=/soft/apps/packages/python-2.7.9/lib ; export LD_LIBRARY_PATH
    else
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/soft/apps/packages/python-2.7.9/lib
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/soft/apps/packages/python-2.7.9/bin:/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/soft/apps/packages/python-2.7.9/bin:/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_12.04-x86_64 case choice
  
  linux-Ubuntu_14.04-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/soft/apps/packages/python-2.7.9/bin:/soft/apps/packages/gcc/gcc-6.2.0/bin:/soft/apps/packages/climate/mpich/3.2/gcc-6.2.0/bin:/soft/apps/packages/hdf5/hdf5-1.8.20/bin:/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/soft/apps/packages/python-2.7.9/bin:/soft/apps/packages/gcc/gcc-6.2.0/bin:/soft/apps/packages/climate/mpich/3.2/gcc-6.2.0/bin:/soft/apps/packages/hdf5/hdf5-1.8.20/bin:/usr/lib/lightdm/lightdm:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append LD_LIBRARY_PATH
    if [ "x${LD_LIBRARY_PATH}" = "x" ]; then
      LD_LIBRARY_PATH=/soft/apps/packages/python-2.7.9/lib:/soft/apps/packages/climate/szip/2.1/gcc-6.2.0/lib/:/soft/apps/packages/gcc/gcc-6.2.0/lib64:/soft/apps/packages/climate/mpich/3.2/gcc-6.2.0/lib:/soft/apps/packages/hdf5/hdf5-1.8.20/lib ; export LD_LIBRARY_PATH
    else
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/soft/apps/packages/python-2.7.9/lib:/soft/apps/packages/climate/szip/2.1/gcc-6.2.0/lib/:/soft/apps/packages/gcc/gcc-6.2.0/lib64:/soft/apps/packages/climate/mpich/3.2/gcc-6.2.0/lib:/soft/apps/packages/hdf5/hdf5-1.8.20/lib
    fi
    
    # translating: set LIBRARY_PATH
    LIBRARY_PATH=/usr/lib/x86_64-linux-gnu ; export LIBRARY_PATH
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/soft/apps/packages/gcc/gcc-6.2.0/share/man:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/soft/apps/packages/gcc/gcc-6.2.0/share/man:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_14.04-x86_64 case choice
  
  linux-Ubuntu_16.04-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_16.04-x86_64 case choice
  
  linux-Ubuntu_6.06.2-ia32)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-Ubuntu_6.06.2-ia32 case choice
  
  linux-Ubuntu_6.06.2-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_6.06.2-x86_64 case choice
  
  linux-Ubuntu_7.04-ia32)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_7.04-ia32 case choice
  
  linux-Ubuntu_7.04-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_7.04-x86_64 case choice
  
  linux-Ubuntu_8.04-ia32)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-Ubuntu_8.04-ia32 case choice
  
  linux-Ubuntu_8.04-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-Ubuntu_8.04-x86_64 case choice
  
  linux-debian_3.0-ia32)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-debian_3.0-ia32 case choice
  
  linux-debian_3.1-ia32)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-debian_3.1-ia32 case choice
  
  linux-debian_4.0-ia32)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-debian_4.0-ia32 case choice
  
  linux-debian_4.0-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-debian_4.0-x86_64 case choice
  
  linux-debian_unstable-x86_64)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/bin:/usr/sbin:/usr/X11R6/bin:/sbin:/usr/local/bin:/usr/local/sbin:/mcs/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/share/man:/usr/local/man:/usr/local/share/man:/usr/X11R6/man:/opt/man:/software/common/adm/man
    fi
  
  ;;
  #  End of linux-debian_unstable-x86_64 case choice
  
  linux-rh73)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-rh73 case choice
  
  linux-tg)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of linux-tg case choice
  
  solaris-2)
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
  
  ;;
  #  End of solaris-2 case choice
  
  solaris-9)
    
    # translating: append LD_LIBRARY_PATH
    if [ "x${LD_LIBRARY_PATH}" = "x" ]; then
      LD_LIBRARY_PATH=/usr/openwin/lib ; export LD_LIBRARY_PATH
    else
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/openwin/lib
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/bin:/usr/openwin/bin:/usr/sbin:/etc:/opt/SUNWspro/bin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/bin:/usr/openwin/bin:/usr/sbin:/etc:/opt/SUNWspro/bin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
    
    # translating: set OPENWINHOME
    OPENWINHOME=/usr/openwin ; export OPENWINHOME
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/usr/openwin/man:/opt/SUNWspro/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/usr/openwin/man:/opt/SUNWspro/man:/software/common/adm/man
    fi
  
  ;;
  #  End of solaris-9 case choice
  
  sun4)
    
    # translating: set OPENWINHOME
    OPENWINHOME=/usr/openwin ; export OPENWINHOME
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/usr/man:/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/usr/man:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/ucb:/bin:/usr/bin:/etc:/usr/etc:/usr/sbin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/ucb:/bin:/usr/bin:/etc:/usr/etc:/usr/sbin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of sun4 case choice
  
  *)
    
    # translating: append MANPATH
    if [ "x${MANPATH}" = "x" ]; then
      MANPATH=/software/common/adm/man ; export MANPATH
    else
      MANPATH=${MANPATH}:/software/common/adm/man
    fi
    
    # translating: append PATH
    if [ "x${PATH}" = "x" ]; then
      PATH=/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin ; export PATH
    else
      PATH=${PATH}:/usr/bin:/usr/sbin:/bin:/sbin:/mcs/bin:/usr/local/bin:/software/common/bin:/soft/apps/bin:/soft/gnu/bin:/soft/com/bin:/soft/adm/bin:${HOME}/bin/${WHATAMI}:${HOME}/bin
    fi
  
  ;;
  #  End of default case choice

esac
# End of case $${ARCH}

# translating: uninitvar PATH
SOFT_REMOVE=/INIT:/bin:/sbin:/usr/local/bin
PATH=`echo $PATH | sed -e "s#:${SOFT_REMOVE}[^:]*##g" -e "s#^${SOFT_REMOVE}[^:]*:\{0,1\}##"`
unset SOFT_REMOVE

# translating: switch ${ARCH}
case ${ARCH} in
  
  irix-6)
  
  ;;
  #  End of irix-6 case choice
  
  solaris-2)
  
  ;;
  #  End of solaris-2 case choice
  
  solaris-9)
  
  ;;
  #  End of solaris-9 case choice

esac
# End of case $${ARCH}

# This is the end of /homes/joneal/.soft.cache.sh
