#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/wk99_power ${npt} /home/ss901165/src/mjo/wk99_power.f90 -lfftpack

