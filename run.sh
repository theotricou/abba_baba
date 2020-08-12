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

  #sed "s/1123581321/$seed/g" $zomp > $dir/zombi_parameters
  cp $zomp $dir/zombi_parameters

  sed "s/1123581321/$seed/g" $simp > $dir/sim_parameters

  singularity run $singularity python /simulation/ZOMBI/Zombi.py T $dir/zombi_parameters $dir

  singularity run $singularity python ${path%run.sh}python3/build_ms_command.py $dir/T/CompleteTree.nwk -p $dir/sim_parameters -s $sampling -o $dir -v

  singularity run $singularity Rscript ${path%run.sh}R/ms_simulation.R $dir

else

  echo
  echo Warning: directory exits, delete before running again or use a different name
  echo

fi

#
# # temp TEST
#
# #!/usr/bin/env Rscript
# # Theo
#
# #shell
#
for i in `seq 1001 1 1001`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  sbatch temp
  rm temp
done
#
# for i in test*;do
#   if [[ `ls $i | wc -l` != 8 ]]; then
#     if  [[ ! `squeue -o %j | grep -P "${i#test}"` ]]; then
#       echo $i
#       rm -rf $i
#     fi
#   fi
# done
#
#
#
#
# rm -rf D*
# for j in sim_*;do
#   cd $j
#   rm -rf D*
#   for i in test*; do
#     sed "1d" $i/data.txt | sed "s/^/$i\t/" >> Dstat
#     sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> D3
#   done
#   cd ..
#   sed "s/^/$j\t/" $j/Dstat >> Dstat
#   sed "s/^/$j\t/" $j/D3 >> D3
# done
#
#
# for j in one_*;do
#   cd $j
#   rm -rf oneD*
#   for i in test*; do
#     sed "1d" $i/data.txt | sed "s/^/$i\t/" >> oneDstat
#     sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> oneD3
#   done
#   cd ..
#   sed "s/^/$j\t/" $j/oneDstat >> oneDstat
#   sed "s/^/$j\t/" $j/oneD3 >> oneD3
# done
#
# sed "s/^/05\t/" oneDstat > oneDstat2
#
#
#
#
#
# # temp R truc
#
# require('ape')
# tree <- read.tree("aletree")
# d<-read.table('ale_transfers')
#
#
# edge<-tree$edge
# spnd<-c(tree$tip.label, tree$node.label)
# wherefrom<-which(spnd==from)
# whereto<-which(spnd==to)
# fromdad<-edge[edge[,2]==wherefrom,1]
# todad<-edge[edge[,2]==whereto,1]
#
#
#
#
# tresh=0.7
# a<-by(d,d[1],function(x) c(sum(x[x[,3]>tresh,3]),nrow(x)))
# b<-as.data.frame(t(do.call(cbind,a)))
# b$names<-as.character(rownames(b))
# rownames(b)<-NULL
# brl<-as.data.frame(cbind(br<-c(tree$tip.label,tree$node.label), c(0,tree$edge.length)))
# brl<-brl[order(brl[,1]),]
# b<-b[order(b$names),]
# b$dist<-brl$V2
# plot(as.numeric(as.character(b$dist)),b$V1, xlim=c(0,0.11))
# plot(as.numeric(as.character(b$dist)),b$V1)
# b[which(b$V1>140),]
# plot(b$V1, b$V2)
