#!/usr/bin/env Rscript
#Theo
# bear data pipeline

#ddl at https://datadryad.org/stash/dataset/doi:10.5061/dryad.cr1496b

shell part

for i in *fa; do
  B="$(cut -d'_' -f2 <<<"$i")"
  mkdir $B
  mv $i $B/
done

for i in */;do
  cd $i
  csplit *.fa '/scaffold/' '{*}'
  cd ..
done

for i in *;do
  cd $i
  for j in xx*; do
    sed -i "s/scaffold.*/$i/g" $j
  done
  cd ..
done


for i in xx*; do
  cat ../Adm1/$i ../235/$i ../Swe/$i ../Uamericanus/$i > ../ali2/$i
done
for i in  xx*; do seaview -convert -sites  $i; done


for i in *fst;do
  sed -i "s/N/-/g" $i
  ~/Downloads/9_softwares/Gblocks_Linux64_0.91b/Gblocks_0.91b/Gblocks $i -b5=n
  snp-sites $i-gb
done


for i in *fst;do
  snp-sites $i
done



for i in *fa; do
  B="$(cut -d'_' -f2 <<<"$i")"
  echo $B
  sed "s/^>/>$B\_/g" $i > $B.fasta
done

rm test
head -n 100 *fasta | sed "s/=.*//g" | sed '/^$/d' >> test
sed "s/N/-/g" test > test2


################################################################################
################################################################################

for i in xx*;do
  seaview -o $i\_n -concatenate $file $i
  rm $file
  msa_split -w 1000000,0 $i\_n --out-root $i 2> log
  file=`tail -n 2 log | grep -ho "xx.*fa"`
done

for i in *fa;do
  snp-sites -c $i
done

################################################################################
################################################################################


library('seqinr')
library(parallel)
files <- list.files(pattern = ".snp_sites.aln")

filtering<-function(x){
  x2=x
  # x2=x[!x %in% c("n", '-')]
  if (length(unique(x2)) != 1){
    if (length(x2)>3){
      if (length(which(table(x2) >= 2)) > 1){
        return(TRUE)
      }else{FALSE}
    }else{FALSE}
  }else{FALSE}
}


for (i in 1:length(files)){
  ali<-as.data.frame(as.matrix(read.alignment(files[i], format = "fasta")))
  print(files[i])
  # fa<-ali[,which(unlist(mclapply(ali,function(x) filtering(x), mc.cores = 6)))]
  fa<-ali[,which(unlist(lapply(ali,function(x) filtering(x))))]
  write.table(fa, paste(files[i], "sites", sep="."), sep = " ", col.names=F,row.names = T, append = F, quote=F)
}
#
#
# a<-as.data.frame(as.matrix(read.alignment('xx99.fst', format = "fasta")))
# dim(a)
# # a2<-a[,which(apply(a,2, function(x) length(unique(x)) != 1))]
# # fa<-a[,which(apply(a,2, function(x) filtering(x)))]
# # dim(fa)
#
# # temp=a[, c(1:10000)]
# # fa<-temp[,which(unlist(mclapply(temp,function(x) filtering(x), mc.cores = 20)))]
# # fa<-temp[,which(unlist(lapply(temp,function(x) filtering(x))))]
# fa<-a[,which(unlist(mclapply(a,function(x) filtering(x), mc.cores = 20)))]
# dim(fa)
#
# write.table(fa, "xx99.fst.sites", sep = " ", col.names=F,row.names = T, append = F, quote=F)
#
#
#
#
#
#
#


Der<-function(x){
  if (x[1] == x[3]){
    return("baba")
  }else if(x[2] == x[3]){
    return("abba")
  }else if(x[2] == x[1]){
    return("bbaa")
  }
}


DD<-function(x){
  d=read.table(x, h=F, row.names=1)
  abba=0
  baba=0
  bbaa=0
  temp=table(unlist(lapply(d, function(x) Der(x))))
  D=((temp['abba']-temp['baba'])/(temp['abba']+temp['baba']))
  return(as.numeric(D))
}

Ds=lapply(files, function(x) DD(x))





Der<-function(x){
  if (x[1] == x[3]){
    return("baba")
  }else if(x[2] == x[3]){
    return("abba")
  }else if(x[2] == x[1]){
    return("bbaa")
  }
}

D_tabler<-function(x){
  temp=table(unlist(lapply(x, function(y) Der(y))))
  return(temp[c("abba", "baba")])
}


files <- list.files(pattern = "_y")
all_data<-lapply(files, function(x) read.table(x, h=F, row.names=1))
all_data[1]

D_table<-lapply(all_data, function(x) D_tabler(x))
truc=do.call('rbind', D_table)


res<-c()
for (i in 1:nrow(truc)){
  truc2<-truc[-i,]
  res[i] = (sum(truc2[,1])-sum(truc2[,2]))/(sum(truc2[,1])+sum(truc2[,2]))
}



The size of blocks should be chosen to be larger than the extent of linkage
disequilibrium. By computing the variance of the D statistic over the sequences
M times leaving each block of the sequence in turn, and then multiplying by M
and taking the square root, we can obtain an approximately normally distributed
standard error using the theory of the jackknife (Reich et al. 2009).

Z = D/sqrt((var(res)*1064))
= 20?
20.35548

D <- 0.103559
n_blocks <- nrow(truc)
D_sd <- sd(res)
D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err



Da=0.079782
Sa=0.009739
Za=8.192120



for (i in lenght(files)){
  file2 = files[-i]
  abba=0
  baba=0
  bbaa=0

}







files <- list.files(pattern = "_y")
for (i in 1:length(files)){
  d=read.table(files[i], h=F, row.names=1)
  print(files[i])
  temp=lapply(d, function(x) Der(x))
}
