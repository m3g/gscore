#!/bin/bash
pdblist=$1    # Name of file containing the model list
logdir=$2     # Name of output directory
nprocs=$3     # Number of processors of your machine
#
# No modifications required after this
#
username=`whoami`
for model in `cat $pdblist`; do
  name=$(basename $model)
  name=${name%.*}
  log=$logdir/$name.log
  wait_for.tcl $username lovoalign $nprocs 0.1
  echo "Running $name"
  lovoalign -p1 $model -pdblist $pdblist -seqfix -skip -nglobal 1 > $log &
done
