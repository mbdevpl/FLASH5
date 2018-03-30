#!/bin/bash

#-------------------------------------------------------------------------------
#
# FlashTest does not communicate success or failure using the error code.
# However, this is how Jenkins determines success or failure of each command
# executed in a Build step.
#
# This script
# 1) sets up a clean test environment,
# 2) runs all desired tests in a single FlashTest run regardless of the success
#    of previous tests,
# 3) checks the root errors file for errors, and
# 4) communicates success/failure via the error code.
#
# For FlashTest to work properly, users of this script must first configure
# compute001 with the correct compilers and libraries using the softenv
# configuration environment.
#
#-------------------------------------------------------------------------------

# Check environment variables
# FLASH_BASE        - full path to root of FLASH directory tree
# FLASHTEST_BASE    - full path to root of FlashTest repository
[ -z "$FLASH_BASE" ]     && { echo "Need to set FLASH_BASE";     exit 1; }
[ -z "$FLASHTEST_BASE" ] && { echo "Need to set FLASHTEST_BASE"; exit 1; }

# Where to write file to (in Jenkins-accessible workspace)
FLASHTEST_OUTPUT=$FLASH_BASE/../TestResults

# Remove old results
RESULTS_DIR=$FLASHTEST_OUTPUT/compute001
if [[ -d  $RESULTS_DIR ]]; then
    rm -rf $RESULTS_DIR
fi

# Build brand-spanking-new version of sfocu
cd $FLASH_BASE/tools/sfocu
make -f Makefile.hello clean
make -f Makefile.hello

# Run all tests in one execution
mkdir -p $FLASHTEST_OUTPUT
$FLASHTEST_BASE/flashTest.py \
                -o $FLASHTEST_OUTPUT \
                -c $FLASH_BASE/sites/compute001.mcs.anl.gov/FlashTest/config \
                -i $FLASH_BASE/sites/compute001.mcs.anl.gov/FlashTest/test.info \
                UnitTest/Grid/AMR/AMReX/2d/Init \
                UnitTest/Grid/AMR/AMReX/2d/Refine \
                UnitTest/Grid/AMR/AMReX/2d/TestCyl \
                Comparison/Sod/PseudoUG/2d/Paramesh/simpleUnsplit \
                Comparison/Sod/AMR/2d/Paramesh/simpleUnsplit \
                Comparison/Sedov/PseudoUG/2d/Paramesh/simpleUnsplit \
                Comparison/Sedov/AMR/2d/Paramesh/simpleUnsplit \
                Comparison/Sod/PseudoUG/2d/Paramesh/unsplit \
                Comparison/Sod/AMR/2d/Paramesh/unsplit

# Confirm output directory contains only one directory, which is the folder
# containing the FlashTest results
if [[ ! -d  $RESULTS_DIR ]]; then
    echo "Expected results directory - $RESULTS_DIR - does not exist"
    exit 2
fi
NDIR=$(ls -d $RESULTS_DIR/* | wc -l)
if [[ $NDIR -ne 1 ]]; then
    echo "Expected only one directory in $RESULTS_DIR"
    echo "Confirm that flashTest.py is only being called once"
    exit 3
fi

# An error occurred if this file is not empty
ERROR_LOG=$(ls -d $RESULTS_DIR/*)/errors
echo
echo "------------------------------------------------------------"
echo "FlashTest Error Log = $ERROR_LOG"
if [[ ! -f $ERROR_LOG ]]; then
    echo "FlashTest error log not found"
    echo "------------------------------------------------------------"
    echo
    exit 4
elif [[ -s $ERROR_LOG ]]; then
    echo "FlashTest reports FAILURE"
    echo
    cat $ERROR_LOG
    echo
    echo "------------------------------------------------------------"
    echo
    exit 5
else
    echo "FlashTest reports SUCCESS"
fi
echo "------------------------------------------------------------"
echo

