
scancel -u tricou
rm -rf test* slurm*
#
for i in `seq 90001 1 90500`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  sbatch temp
  #bash temp
  rm temp
done

grep -P 'N2\tP' test[34][0-9]/data.txt | column -t | grep 'N2'



rm data_D.txt data_Dfoil.txt
files=(test*)
# get the 1st element of the array
first=${files[0]}
# remove the 1st element of from array
files=("${files[@]:1}")
# process the 1st element
sed "s/^/`echo $first`\t/" $first/data.txt > data_D.txt
sed "s/^/`echo $first`\t/" $first/data_Dfoil.txt > data_Dfoil.txt
name=`echo $first`
sed -i -e "0,/$name/ s/$name/test\t/" data_D.txt
sed -i -e "0,/$name/ s/$name/test\t/" data_Dfoil.txt
# process the remaining elements
for i in "${files[@]}"; do
  sed "1d" $i/data.txt | sed "s/^/$i\t/" >> data_D.txt
  sed "1d" $i/data_Dfoil.txt | sed "s/^/$i\t/" >> data_Dfoil.txt
done
head data_D.txt
head data_Dfoil.txt

grep "Error" data_Dfoil.txt |sort |  uniq
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




for i in test*;do
  if [[ `ls $i | wc -l` != 7 ]]; then
    if  [[ ! `squeue -o %j | grep -P "${i#test}"` ]]; then
      echo $i
      rm -rf $i
    fi
  fi
done
for i in `seq 90001 1 90500`; do
  if [ ! -d test$i ]; then
    sed "s/aaaa/$i/g" run_slurm.sh > temp
    sbatch temp
    # bash temp
    rm temp
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

for i in `seq 10001 1 10500`; do
  if [ ! -d test$i ]; then
    sed "s/aaaa/$i/g" run_slurm.sh > temp
    sbatch temp
    # bash temp
    rm temp
  fi
done


for i in test*; do
  if (( `ls $i | wc -l` != 7));then
    echo $i
    rm -rf $i
    sed "s/aaaa/${i:4}/g" run_slurm.sh > temp
    sbatch temp
    rm temp
  fi
done
