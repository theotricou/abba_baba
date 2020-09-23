
scancel -u tricou
rm -rf test* slurm*
#
for i in `seq 1001 1 1001`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  # sbatch temp
  # rm temp
done


rm data_D.txt data_D3.txt
files=(test*)
# get the 1st element of the array
first=${files[0]}
# remove the 1st element of from array
files=("${files[@]:1}")
# process the 1st element
sed "s/^/`echo $first`\t/" $first/data.txt > data_D.txt
sed "s/^/`echo $first`\t/" $first/data_D3.txt > data_D3.txt
name=`echo $first`
sed -i -e "0,/$name/ s/$name/test\t/" data_D.txt
sed -i -e "0,/$name/ s/$name/test\t/" data_D3.txt
# process the remaining elements
for i in "${files[@]}"; do
  sed "1d" $i/data.txt | sed "s/^/$i\t/" >> data_D.txt
  sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> data_D3.txt
done
head data_D.txt
head data_D3.txt


# colles + proba

for i in tes*; do
  if [ ! -f $i/col_pro.txt ]; then
    sed "s/aaaa/$i/g" /beegfs/data/tricou/abba_baba/coller_slurm.sh > temp
    sbatch temp
    rm temp
  fi
done



for i in test*; do
  echo $i `grep -oh " at .*" $i/ms_command.R | awk '{print $2}'`
done > TR_time




# for i in test*;do
  if [[ `ls $i | wc -l` != 8 ]]; then
    if  [[ ! `squeue -o %j | grep -P "${i#test}"` ]]; then
      echo $i
      rm -rf $i
    fi
  fi
done




rm -rf D*
for j in sim_*;do
  cd $j
  rm -rf D*
  for i in test*; do
    sed "1d" $i/data.txt | sed "s/^/$i\t/" >> Dstat
    sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> D3
  done
  cd ..
  sed "s/^/$j\t/" $j/Dstat >> Dstat
  sed "s/^/$j\t/" $j/D3 >> D3
done


for j in one_*;do
  cd $j
  rm -rf oneD*
  for i in test*; do
    sed "1d" $i/data.txt | sed "s/^/$i\t/" >> oneDstat
    sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> oneD3
  done
  cd ..
  sed "s/^/$j\t/" $j/oneDstat >> oneDstat
  sed "s/^/$j\t/" $j/oneD3 >> oneD3
done

sed "s/^/05\t/" oneDstat > oneDstat2


for i in test*; do
  cd $i
  sbatch /beegfs/data/tricou/abba_baba/sisterer.sh
  cd ..
done

for i in test*;do
  if [ ! -f $i/data_s.txt ]; then
    cd $i
    sbatch /beegfs/data/tricou/abba_baba/sisterer.sh
    cd ..
  fi
done



# temp R truc

require('ape')
tree <- read.tree("aletree")
d<-read.table('ale_transfers')


edge<-tree$edge
spnd<-c(tree$tip.label, tree$node.label)
wherefrom<-which(spnd==from)
whereto<-which(spnd==to)
fromdad<-edge[edge[,2]==wherefrom,1]
todad<-edge[edge[,2]==whereto,1]




tresh=0.7
a<-by(d,d[1],function(x) c(sum(x[x[,3]>tresh,3]),nrow(x)))
b<-as.data.frame(t(do.call(cbind,a)))
b$names<-as.character(rownames(b))
rownames(b)<-NULL
brl<-as.data.frame(cbind(br<-c(tree$tip.label,tree$node.label), c(0,tree$edge.length)))
brl<-brl[order(brl[,1]),]
b<-b[order(b$names),]
b$dist<-brl$V2
plot(as.numeric(as.character(b$dist)),b$V1, xlim=c(0,0.11))
plot(as.numeric(as.character(b$dist)),b$V1)
b[which(b$V1>140),]
plot(b$V1, b$V2)
