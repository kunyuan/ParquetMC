#!/bin/bash
#set -x
########################################################################
# SUN Grid Engine job wrapper
########################################################################
#$ -N renorm
#$ -q i08m3
#$ -j y
#$ -M chenkun0228@gmail.com
#$ -m e
#$ -v WIEN_DMFT_ROOT,WIENROOT,LD_LIBRARY_PATH,PATH
########################################################################
# DON'T remove the following line!
source $TMPDIR/sge_init.sh
########################################################################
export SMPD_OPTION_NO_DYNAMIC_HOSTS=1
export OMP_NUM_THREADS=1
export PATH=.:$PATH
export MODULEPATH=/opt/apps/modulefiles:/opt/intel/modulefiles:/opt/pgi/modulefiles:/opt/gnu/modulefiles:/opt/sw/modulefiles

# Adjust modules if needed
module load intel/ompi

for var in {0..31}
do
    ./feyncalc.exe < ./_in$var >& ./_out$var
done
