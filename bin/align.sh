#!/bin/bash
pdblist=$1    # Name of file containing the model list
logdir=$2     # Name of output directory
nprocs=$3     # Number of processors of your machine
#
# No modifications required after this
#
username=`whoami`
for file in `cat $pdblist`; do
  name=$(basename $model)
  name=${name%.*}
  log=$logdir/$name.log
  wait_for.tcl $whoami lovoalign $nprocs
  lovoalign -pdblist $pdblist -seqfix -skip > $log &
done
