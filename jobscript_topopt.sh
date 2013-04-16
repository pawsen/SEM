#!/bin/sh

# run with
# qsub jobscript.sh
# http://beige.ucs.indiana.edu/I590/node35.html
# qstat: R: running, C:completed, Q: queued, S: suspended,
#        H:held this means that it is not going to run until it is released
# showq:alle jobs

# 1 nodes x 12 cores
#PBS -l nodes=4:ppn=12
#PBS -q topopt

# embedded options to qsub - start with #PBS
# -- our name ---
#PBS -N SEM

#PBS -o stout/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e stout/$PBS_JOBNAME.$PBS_JOBID.err
# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

# does not have any effect on MPI, but if you mix mpi with openmp it will have
# an effect
# export OMP_NUM_THREADS=8

source /appl/htopopt/RH61/set_env.sh
MPIRUN=/appl/htopopt/RH61/openmpi-1.4.5/bin/mpirun
# $MPIRUN --bind-to-socket GOL_P

# $MPIRUN -np 8 -mca btl tcp,self GOL_P
# Set enviroment
module load VTK

arr=(2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48)
for nodes in ${arr[@]}
do
    $MPIRUN -np $nodes -mca btl tcp,self SEM_p
    #echo $nodes
done
