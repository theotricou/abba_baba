#!/bin/bash
# Theo
#run complete simulation for ABBA BABA test.
path=$0
singularity=$1
dir=$2
seed=$3
zomp=$4
simp=$5
sampling=$6

if [ ! -d "$dir" ]; then
  mkdir $dir

  # switch command for single tree or random trees simulation
  sed "s/1123581321/$seed/g" $zomp > $dir/zombi_parameters
  # cp $zomp $dir/zombi_parameters
  ##

  sed "s/1123581321/$seed/g" $simp > $dir/sim_parameters

  singularity run $singularity python /simulation/ZOMBI/Zombi.py T $dir/zombi_parameters $dir

  singularity run $singularity python ${path%run.sh}python3/build_ms_command_D3.py $dir/T/CompleteTree.nwk -p $dir/sim_parameters -s $sampling -o $dir -v

  singularity run $singularity Rscript ${path%run.sh}R/ms_simulation_D3.R $dir

else

  echo
  echo Warning: directory exits, delete before running again or use a different name
  echo

fi
