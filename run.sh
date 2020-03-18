#!/bin/bash
# Theo
#run complete simulation for ABBA BABA test.
path=$0
singularity=$1
dir=$2
seed=$3
zomp=$4
simp=$5

if [ ! -d "$dir" ]; then
  mkdir $dir

  sed "s/1123581321/$seed/g" $zomp > $dir/zombi_parameters
  sed "s/1123581321/$seed/g" $simp > $dir/sim_parameters

  singularity run $singularity python /simulation/ZOMBI/Zombi.py T $dir/zombi_parameters $dir

  singularity run $singularity python ${path%run.sh}python3/build_ms_command.py $dir/T/CompleteTree.nwk -p $dir/sim_parameters -o $dir -v

  singularity run $singularity Rscript ${path%run.sh}R/ms_simulation.R $dir

else

  echo
  echo Warning: directory exits, delete before running again or use a different name
  echo

fi



for j in div*; do
  cd $j
  for i in `seq 666 1 765`; do
    sed "s/aaaa/$i/g" run_slurm.sh > temp
    sbatch temp
    rm temp
  done
  cd ..
done

rm second
for i in tes*/data*; do
  sed "1d" $i >> second
done


# test on N2 donor
grep " [0-9]*@[1-9]" tes*/ms_command.R


for i in `seq 1000 1 1000`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  # sbatch temp
  # rm temp
done


for i in `seq 1300 1 1499`;do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  sbatch temp
  rm temp
done



for i in `seq 1100 1 1299`;do
  if [ ! -f test$i/ms_command.R ];then
    rm -rf test$i
    sed "s/aaaa/$i/g" run_slurm.sh > temp
    sbatch temp
    rm temp
  fi
done
