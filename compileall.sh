#!/usr/bin/ksh
#
# Note, if ksh is not available, you can use
# another shell like zsh
#
##################################################
##################################################
# This script will compile all executable programs
# that are required by the FMM based tomographic
# inversion routine FMTT. Additional programs for
# plotting, computing model roughness etc. are
# also compiled. See accompanying documentation
# for further details
##################################################
##################################################
#
# YOU MUST SPECIFY A FORTRAN 90 COMPILER!!!!
# Compiler options can be included at this stage
# if deemed necessary by the user.
#
##################################################
F90="gfortran -O3"
##################################################
# NOTE: Compilation of the code "aktsurf", which
# computes ak135 traveltimes to the base of the
# local 3-D model (in conjunction with "itimes"),
# is also carried out by this script. However, you
# must ensure that the Makefile in the subdirectory
# "source/aktimes" points to an appropriate FORTRAN
# compiler. Ideally, use the same compiler that
# you use above.
##################################################
#
# Enter the directory called "source" and compile
# all programs.
#
##################################################
cd source_mod
$F90 -o fm3dt fm3dt.f90 -g
echo " "
echo "Compilation of fm3dt complete"
$F90 -o gmtslicet gmtslicet.f90 -g
echo " "
echo "Compilation of gmtslicet complete"
$F90 -o grid3dtg grid3dtg.f90 -g
echo " "
echo "Compilation of grid3dtg complete"
$F90 -o itimes itimes.f90 -g
echo " "
echo "Compilation of itimes complete"
$F90 -o misfitt misfitt.f90 -g
echo " "
echo "Compilation of misfitt complete"
$F90 -o residualst residualst.f90 -g
echo " "
echo "Compilation of residualst complete"
$F90 -o resplott resplott.f90 -g
echo " "
echo "Compilation of resplott complete"
$F90 -o subinv subinv.f90 -g
echo " "
echo "Compilation of subinv complete"
$F90 -o syntht syntht.f90 -g
echo " "
echo "Compilation of syntht complete"
$F90 -g -c lsmr.f90 -o lsmr.o -O3
$F90 -g -c fmttinv.f90 -o fmttinv.o -O3
$F90 fmttinv.o lsmr.o -o fmttinv -O3
##################################################
#
# Move all executables to directory ../bin
#
##################################################
rm *.o *.mod
mv fm3dt gmtslicet grid3dtg itimes misfitt fmttinv ../bin
mv residualst resplott subinv syntht ../bin
##################################################
#
# Enter subdirectory "aktimes" and compile program
# for generating traveltime tables
#
##################################################
cd aktimes
make all
make clean
./remodl
./setbrn
##################################################
#
# Move "aktsurf" to ../../bin and copy binary
# traveltime tables to ../inputfiles
#
##################################################
mv aktsurf ../../bin
cp ak135.tbl ak135.hed ../inputfiles
##################################################
cd ../..
echo " "
echo "Compilation complete"
