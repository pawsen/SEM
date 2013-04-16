#!/bin/sh

# run with
# qsub jobscript.sh
# http://beige.ucs.indiana.edu/I590/node35.html
# qstat: R: running, C:completed, Q: queued, S: suspended,
#        H:held this means that it is not going to run until it is released
# showq:alle jobs

# embedded options to qsub - start with #PBS
# -- our name ---
#PBS -N SEM

#PBS -o stout/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e stout/$PBS_JOBNAME.$PBS_JOBID.err
# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

# 4 nodes x 8 cores
#PBS -l nodes=1:ppn=1
#PBS -q hpc

# does not have any effect on MPI, but if you mix mpi with openmp it will have
# an effect
# export OMP_NUM_THREADS=8

source /appl/htopopt/RH61/set_env.sh
MPIRUN=/appl/htopopt/RH61/openmpi-1.4.5/bin/mpirun
# $MPIRUN --bind-to-socket GOL_P

# $MPIRUN -np 8 -mca btl tcp,self GOL_P
# Set enviroment
module load VTK

./SEM_gbar

# arr=(2 4 8 16 32)
# for nodes in ${arr[@]}
# do
#     $MPIRUN -np $nodes -mca btl tcp,self GOL_P
#     #echo $nodes
# done
