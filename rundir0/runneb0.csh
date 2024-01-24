#!/bin/csh -f
#PBS -l select=1:ncpus=2:mpiprocs=1:ompthreads=1:jobtype=core
#PBS -l walltime=72:00:00 

module purge

if ($?PBS_O_WORKDIR) then
cd ${PBS_O_WORKDIR}
endif

source /local/apl/lx/g16a03/g16/bsd/g16.login

mkdir /work/users/to1/tmp.$$
set WORK=/work/users/to1/tmp.$$

python3 /home/users/to1/NEBprog/neb_starter.py neb.com ${WORK}

exit 0

