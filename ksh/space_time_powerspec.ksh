#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/space_time_powerspec ${npt} /home/ss901165/src/mjo/space_time_powerspec.f90 -lfftpack
