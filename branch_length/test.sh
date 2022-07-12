rm *

echo SEP > AAA
for i in {1..1}
do
  /home/theo/Downloads/9_softwares/ms.folder/msdir/ms \
    4 1000 -I 4 1 1 1 1 -t 10 -ej 2 3 2 -ej 3 2 1 -ej 7 4 1 -T |
    tail -n +4 | grep -v // > a$i.txt
  # Simulate sequences from trees using Seq-Gen.
  /home/theo/Downloads/9_softwares/Seq-Gen-1.3.4/source/seq-gen \
   -m HKY -l 1000 -s 0.02 -of -x AAA < a$i.txt > b$i.txt
done
cat a* | grep ",(3:......("
csplit --suppress-matched b1.txt '/^$/' {*}
rm xx100
for i in xx*; do
  sed -i "s/SEP//" $i
done
rsync -e ssh -avz xx* tricou@pbil-deb:/beegfs/data/tricou/abba_baba/temp/

sss a1.txt


library("seqinr")

sites<-function(x){
  if (length(unique(c(x))) == 2){
    if ( (x[1] == x[2]) & (x[3] == x[4]) ){"bbaa"}
    else if ( (x[1] == x[2]) & (x[1] == x[3]) ){"baaa"}
    else if ( (x[1] == x[2]) & (x[1] == x[4]) ){"abaa"}
    else if ( (x[1] == x[3]) & (x[1] == x[4]) ){"aaba"}
  }
}

counter<-function(fasta_file){
  file<-read.fasta(fasta_file)
  second<-attr(file[2], "name")
  if (second == "1"){To = "ABC"}
  if (second == "3"){To = "BCA"}
  if (second == "2"){To = "ACB"}
  fasta<-as.data.frame(lapply(file, function(x) x[1:1000]))
  # apply(fasta, 1, function(x) sites(x))
  snp<-table(unlist(apply(fasta, 1, function(x) sites(x))))
  TT1<-as.numeric((((snp['baaa'] + snp['abaa']) / 2) + snp['bbaa']) / 1000)
  TT2<-as.numeric(((snp['baaa'] + snp['abaa']) / 2) / 1000)
  return(c("Topology" = To, "TT1" = TT1, "TT2" = TT2))
}
std <- function(x) sd(x)/sqrt(length(x))

files=Sys.glob("xx*")

temp <- lapply(files, counter)
RES<-as.data.frame(do.call('rbind',temp))
RES$TT1<-as.numeric(as.character(RES$TT1))
RES$TT2<-as.numeric(as.character(RES$TT2))
par(mfrow=c(1,2))
boxplot(RES$TT1 ~ RES$Topology)
boxplot(RES$TT2 ~ RES$Topology)


res<-as.data.frame(do.call('rbind',by(RES, RES$Topology, function(x) {
  res=c("mT1"=mean(x$TT1), "mT2"=mean(x$TT2), "seT1"=std(x$TT1), "seT2"=std(x$TT2))
  return(res)
})))


plot(res$mT1, ylim=c(min(res$mT2), max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n',
  ylab ="", xlab= "",
  main= "all migrations")
axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
points(res$mT2, pch = 19 , cex = 2 , col='blue')
legend("bottomright", legend = c("T1", "T2"), col = c("red", "blue"),
pch = c(19,19), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black",
horiz = F , inset = c(0.1, 0.1))




rm(list=ls())

# library("Biostrings")
# library(parallel)



TT_count<-function(file){
  fastaFile <- readDNAStringSet(file)
  if (all(names(fastaFile) == c("1","2","3","4"))) {To="ABC"}
  if (all(names(fastaFile) == c("1","3","2","4"))) {To="ACB"}
  if (all(names(fastaFile) == c("1","4","2","3"))) {To="BCA"}
  pa=unlist(lapply(seq(1,1000), function(x) {
    site=c(as.character(fastaFile[[1]][x]),as.character(fastaFile[[2]][x]),
           as.character(fastaFile[[3]][x]),as.character(fastaFile[[4]][x]))
    if (length(unique(site)) == 2) {
      if (2 %in% table(site)){
        if (site[1] == site[2]){"bbaa"}
        else{"other"}
      }else if (1 %in% table(site)){
        temp = table(site)
        num=which(site == names(which(temp == 1)))
        if (num == 2){"aaba"}
        else if (num == 3){"abaa"}
        else if (num == 4){"baaa"}
        else{"other"}
      }else{"other"}
    }else{"other"}
  }))
  T2=(1 / 1000) *  ((length(which(pa == 'abaa')) + length(which(pa == 'baaa')))/2)
  T1=(1 / 1000) * (((length(which(pa == 'abaa')) + length(which(pa == 'baaa')))/2) + length(which(pa == 'bbaa')))
  res=c("Topology"=To, "T1" = T1, "T2" = T2,
    "baaa" = length(which(pa == 'baaa')),
    "abaa" = length(which(pa == 'abaa')),
    "aaba" = length(which(pa == 'aaba')),
    "bbaa" = length(which(pa == 'bbaa')),
    "file"=file)
  return(res)
}
std <- function(x) sd(x)/sqrt(length(x))

file =files[6]

files=Sys.glob("xx*")
# lapply(files, function(x) TT_count(x))

TT<-as.data.frame(do.call('rbind',mclapply(files, function(x) TT_count(x), mc.cores = 8)))
TT$T1<-as.numeric(as.character(TT$T1))
TT$T2<-as.numeric(as.character(TT$T2))
res<-as.data.frame(do.call('rbind',by(TT, TT$Topology, function(x) {
  res=c("mT1"=mean(x$T1), "mT2"=mean(x$T2), "seT1"=std(x$T1), "seT2"=std(x$T2))
  return(res)
})))
par(mfrow=c(1,3))
boxplot(TT$T1 ~ TT$Topology)
boxplot(TT$T2 ~ TT$Topology)
plot(res$mT1, ylim=c(min(res$mT2),max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n')
axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
points(res$mT2, pch = 19 , cex = 2 , col='blue')
table(TT$Topology)
legend("bottomright",legend = c("T1", "T2"),col = c("red", "blue"),
  pch = c(19,19), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black",
  horiz = F , inset = c(0.1, 0.1))

TT[which(TT$Topology == "ACB"), ]


fastaFile <- readDNAStringSet(files[83])
fastaFile
TT_count(files[83])

fastaFile <- readDNAStringSet(files[70])
fastaFile
TT_count(files[70])

pa=unlist(lapply(seq(1,1000), function(x) {
  site=c(as.character(fastaFile[[1]][x]),as.character(fastaFile[[2]][x]),
         as.character(fastaFile[[3]][x]),as.character(fastaFile[[4]][x]))
  if (length(unique(site)) == 2) {
    if (2 %in% table(site)){
      if (site[1] == site[2]){"bbaa"}
      else{"other"}
    }else if (1 %in% table(site)){
      temp=table(site)
      num=which(site == names(which(temp == 1)))
      if (num == 2){"aaba"}
      else if (num == 3){"abaa"}
      else if (num == 4){"baaa"}
      else{"other"}
    }else{"other"}
  }else{"other"}
}))

as.data.frame(do.call('rbind',lapply(seq(1,1000), function(x) {
  site=c(as.character(fastaFile[[1]][x]),as.character(fastaFile[[2]][x]),
         as.character(fastaFile[[3]][x]),as.character(fastaFile[[4]][x]))
  if (length(unique(site)) == 2) {
    if (2 %in% table(site)){
      if (site[1] == site[2]) {site}
    }}
})))


as.data.frame(do.call('rbind',lapply(seq(1,1000), function(x) {
  site=c(as.character(fastaFile[[1]][x]),as.character(fastaFile[[2]][x]),
         as.character(fastaFile[[3]][x]),as.character(fastaFile[[4]][x]))
  if (length(unique(site)) == 2) {
    if (1 %in% table(site)){
      temp=table(site)
      num=which(site == names(which(temp == 1)))
if (num == 4){site}
}
}
})))
