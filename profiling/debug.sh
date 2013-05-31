#!/bin/bash
#
# debug.sh - a PBS script to start a debugging session with TotalView
#
# submit as:  COMP=/comp_env PROG=./your_prog qsub [options] debug.sh
#
# where comp_env is your MPI compiler environment, e.g. gcc, studio, or intel
#
# written by: Bernd Dammann <support@hpc.dtu.dk> , 02/2012
#
# name of the job
# #PBS -N TVdebug
# we need to forward the DISPLAY, and the name of our program to debug
#PBS -v DISPLAY,PROG,COMP
# keep the log files local
#PBS  -k oe
# we should be done within one hour
#PBS -l walltime=1:00:00
# start MPI procs on 2 nodes, 8 procs on each
#PBS -l nodes=2:ppn=8
#
# change into work directory
cd $PBS_O_WORKDIR


# load the right MPI support module
#
module load mpi${COMP}

# start TotalView, with $PROG, using OpenMPI, on PBS_NP processors
#
totalview $PROG -mpi "Open MPI" -np $PBS_NP -no_show_startup_parameters 

