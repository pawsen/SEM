#!/bin/sh

# run with
# qsub jobscript.sh
# http://beige.ucs.indiana.edu/I590/node35.html
# qstat: R: running, C:completed, Q: queued, S: suspended,
#        H:held this means that it is not going to run until it is released
# showq:alle jobs

# 1 nodes x 12 cores
#PBS -l nodes=7:ppn=12
#PBS -q topopt
#PBS -l walltime=12:00:00

# embedded options to qsub - start with #PBS
# -- our name ---
#PBS -N SEM

#PBS -o stout/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e stout/$PBS_JOBNAME.$PBS_JOBID.err
# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodefile
sort -u nodefile > mynodes

# does not have any effect on MPI, but if you mix mpi with openmp it will have
# an effect
# export OMP_NUM_THREADS=8

source /appl/htopopt/RH61/set_env.sh
MPIRUN=/appl/htopopt/RH61/openmpi-1.4.5/bin/mpirun

# Set enviroment
module load VTK

JID="paraview"

arr=( 6  12  24  36  48)
arr=(81)
for nodes in ${arr[@]}
do
    $MPIRUN -np $nodes -machinefile ./mynodes --bysocket --bind-to-socket --report-bindings ./SEM_p _${JID}
    #echo $nodes
done
