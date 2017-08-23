#! /bin/bash -l

## LoadLeveler script to submit 2 node, 4 task MPI program ...

# @ job_type = MPICH
# @ class = LONG
# @ account_no = NONE
# @ node = 2
# @ tasks_per_node = 3
# @ wall_clock_limit = 1:00:00
# @ notification = always
# @ notify_user = vbrendel@indiana.edu
# @ environment = COPY_ALL;
# @ output = gsq.$(cluster).$(process).out
# @ error = gsq.$(cluster).$(process).err
# @ queue

## Users should always cd into their execution directory due to
## a bug within LoadLeveler in dealing with the initialdir keyword.

cd ${HOME}/src/GENESEQER/mpidata


## Use mpirun to execute your MPICH program in parallel;
## $LOADL_TOTAL_TASKS and $LOADL_HOSTFILE are defined by
## LoadLeveler for jobs of type MPICH.

#
# Put your Job commands here.
#
#------------------------------------------------
# Specify binary and arguments here:
#
#
#CHANGE THIS (location of binary, arguments to GeneSeqerMPI):

PROG="../bin/GeneSeqerMPI"
PROGARGS="-s Arabidopsis -d ../data/ATest1 ../data/ATest2 -O ./out.gbA.MPI6 -g ../data/U89959"

mpirun -np $LOADL_TOTAL_TASKS -machinefile $LOADL_HOSTFILE $PROG $PROGARGS

#------------------------------------------------
