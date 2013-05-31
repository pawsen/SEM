#!/bin/bash
#
# collect.sh - a PBS script to start a collect session (Studio # analyzer)
#              for MPI programs (multi-node)
#
# submit as:  PROG=./your_prog qsub [options] collect.sh
#
# written by: Bernd Dammann <support@hpc.dtu.dk> , 02/2012
#
# name of the job
#PBS -N collect
# we need to forward the name of our program to debug
#PBS -v PROG
# keep the log files local
#PBS  -k oe
# we should be done within one hour
#PBS -l walltime=1:00:00
# the number of cores we need
#PBS -l procs=8
#
# change into work directory
cd $PBS_O_WORKDIR


# load the right compiler and MPI support module
#
module load studio
module load mpi/studio

PNAME=`echo $PROG | sed 's/ .*$//'` 
STR=`basename $PNAME`
JID=`echo $PBS_JOBID | sed 's/\..*$//'`
# start collect with $PROG, using OpenMPI (aka CT8.1)
#
collect -M OMPT -o collect_${STR}_${JID}.er mpirun -- $PROG 

