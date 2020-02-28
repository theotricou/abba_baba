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

squeue  | grep launch | awk '{print $1}' > ll

while read -r line; do
  scontrol release $line
done < ll



for i in `seq 10000 1 10099`; do
  sed "s/SEED/$i/g" ../rep.sh > temp$i
  sbatch temp$i && rm temp$i
done



dir=666
seaview test\_$dir/spe_tree &
more slurm.sim\_$dir*out

R
a = read.table('test_666/data.txt', h =T)
b =  a[a$Pvalue < 0.001,]

nrow(b)
nrow(b[b$type_D == "P3" & b$type_R == "P1",])
nrow(b[b$type_D == "P3" & b$type_R == "P2",])
nrow(b[b$type_D == "P1" & b$type_R == "P3",])
nrow(b[b$type_D == "P2" & b$type_R == "P3",])
nrow(b[b$type_D == "P4",])
nrow(b[b$type_R == "P4",])


nrow(b[b$type_D == "N3",])
nrow(b[b$type_D == "N2",])
nrow(b[b$type_D == "N1",])
