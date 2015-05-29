#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/wk99_smoothing ${npt} /home/ss901165/src/mjo/wk99_smoothing.f90 -lfftpack

