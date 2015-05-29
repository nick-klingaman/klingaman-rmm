#!/bin/ksh

f90=f90n
npt=${HOME}/src/nicks_subroutines.f90

${f90} -o /home/ss901165/src/mjo/create_segments ${npt} /home/ss901165/src/mjo/create_segments.f90 -lfftpack

