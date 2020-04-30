#!/usr/bin/env Rscript
# Theo

d<-read.table('data.txt',h=T)
d$Donor = as.character(d$Donor)
d$Recip = as.character(d$Recip)
# positive_eve<-d[d$Donor == c('P3') & d$Recip == 'P2'|
#   d$Donor == c('P3') & d$Recip == 'P1'|
#   d$Recip == c('P3') & d$Donor == 'P1'|
#   d$Recip == c('P3') & d$Donor == 'P2',]
signi_D<-d[d$Pvalue<=0.05,]
signi_H0<-signi_D[(signi_D$Donor == c('P3') & signi_D$Recip == 'P1') |
  (signi_D$Donor == c('P3') & signi_D$Recip == 'P2') |
  (signi_D$Recip == c('P3') & signi_D$Donor == 'P1') |
  (signi_D$Recip == c('P3') & signi_D$Donor == 'P2'),]
posD=nrow(signi_H0)/nrow(signi_D)*100
errD=100-posD


d3<-read.table('data_D3.txt',h=T)
d3$Donor = as.character(d3$Donor)
d3$Recip = as.character(d3$Recip)
plot(d3$D3)

d3[which(d3$D3 == min(d3$D3)),]
d3[which(d3$D3 == max(d3$D3)),]

df <-d3[order(d3$Donor, d3$Recip, d3$D3),]
plot(df$D3,pch=19,cex=0.5, col="red")


expected_non_nul3<-which(!d3$Recip == "O" &
!d3$Recip == "N1"&
!(d3$Donor == "P1" & d3$Recip== "P2")&
!(d3$Donor == "P2" & d3$Recip== "P1"))

d3_non_nul<-d3[expected_non_nul3,]
d3_nul<-d3[setdiff(seq(1,nrow(d3)), expected_non_nul3),]

ks.test(d3_non_nul[1,]$D3, d3_nul$D3)


library("fdrtool")
ss = seq(1,length(d3_nul$D3))
W1 <- lm(d3_nul$D3 ~ ss)
CR <- (d3$D3-W1$fitted.values)/summary(W1)$sigma
# mydata$V6 <- (mydata$V3-W1$fitted.values)/summary(W1)$sigma
normCR <- 1-pnorm(CR)
adjustnormCR <- p.adjust(normCR, "fdr")
d3[as.numeric(as.character(names(adjustnormCR[adjustnormCR <= 0.05]))),]


nrow(d3_non_nul)
nrow(d3_nul)

positive_eve3<-d3[d3$Donor == c('P3') & d3$Recip == 'P1'|
  d3$Donor == c('P3') & d3$Recip == 'P2'|
  d3$Recip == c('P3') & d3$Donor == 'P2'|
  d3$Recip == c('P3') & d3$Donor == 'P1',]

d3_alt<-d3[(d3$Donor == "O" & d3$Recip== "P2")|
  (d3$Donor == "O" & d3$Recip== "P1"),]




expected_nul3<-d3[d3$Recip==c("O","N1")|
d3$Donor == "P1" & d3$Recip== "P2"|
d3$Donor == "P2" & d3$Recip== "P1",]

nrow(expected_nul3)


expected_non_nul3<-
d3[!d3$Recip == "O" &
!d3$Recip == "N1"&
!d3$Donor == "N1"&
!(d3$Donor == "P1" & d3$Recip== "P2")&
!(d3$Donor == "P2" & d3$Recip== "P1"),]


nrow(expected_non_nul3)

abline(h = quantile(expected_nul3$D3, probs = c(0.025, 0.975)), col = 'red')

  library("fdrtool")

  W1 <- lm(V3~V2, data=mydata)
  #### W1
  ##centré réduit
  CR <- (mydata$V3-W1$fitted.values)/summary(W1)$sigma
  mydata$V6 <- (mydata$V3-W1$fitted.values)/summary(W1)$sigma
  mydata$V7 <- 1-pnorm(mydata$V6)
  mydata$V8 <- p.adjust(mydata$V7, "fdr")
  mydata[which(mydata$V5 == br),8]





expected_nul3<-d3[
  d3$Recip=="O" &
  d3$Recip=="N1"
  d3$Donor == 'P2'&d3$Recip == 'P1'|
  d3$Donor == 'P1'&d3$Recip == 'P2',]


expected_non_nul3<-d3[!d3$Recip=="" & !d3$Recip=="N1",]


positive_eve3<-d3[d3$Donor == c('P3') & d3$Recip == 'P1'|
  d3$Donor == c('P3') & d3$Recip == 'P2'|
  d3$Recip == c('P3') & d3$Donor == 'P2'|
  d3$Recip == c('P3') & d3$Donor == 'P1',]

signi_D3<-d3[d3$Pvalue<=0.5,]
posD3=nrow(positive_eve3)/nrow(signi_D3)*100
errD3=100-posD3




################################################################################
################################################################################
################################################################################

for i in `seq 1501 1 1700`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  sbatch temp
  rm temp
done

rm -rf test15* test16* test1700 D*
for i in test*; do
  sed "1d" $i/data.txt | sed "s/^/$i\t/" >> Dstat
  sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> D3
done

R


# ev_base=sum(P2P3,P3P2,P1P3,P3P1)
# ev_alt=sum(N2P2,N2P1,OP1,OP2)
#
# err = ev_alt / (ev_alt+ev_base)

expected_error_rate<-function(temp){
  n_H0<-nrow(temp[temp$V13 == "P3P2" |
  temp$V13 == "P2P3"|
  temp$V13 == "P3P1"|
  temp$V13 == "P1P3",])
  n_alt<-nrow(temp[temp$V13 == "N2P2" |
  temp$V13 == "OP1"|
  temp$V13 == "N2P1"|
  temp$V13 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)

  temp$V14<-(p.adjust(temp$V9,"fdr"))
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V14<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V13 == "P3P2" |
  signi_D$V13 == "P2P3"|
  signi_D$V13 == "P3P1"|
  signi_D$V13 == "P1P3",])
  err_all<- 1 - (true_pos/nrow(signi_D))

  alt=nrow(signi_D[signi_D$V13 == "OP2" |
  signi_D$V13 == "OP1"|
  signi_D$V13 == "N2P1"|
  signi_D$V13 == "N2P2",])
  err_alt=alt/nrow(signi_D)
  data=c("expected_error" = exp_err,"observed_error" = err_alt)
  return(data)
}

ds00=read.table('Dstat_00',h=F)
ds01=read.table('Dstat_01',h=F)
ds05=read.table('Dstat_05',h=F)
ds09=read.table('Dstat_09',h=F)
ds00$V13<-as.factor(paste(ds00$V11,ds00$V12, sep=""))
ds01$V13<-as.factor(paste(ds01$V11,ds01$V12, sep=""))
ds05$V13<-as.factor(paste(ds05$V11,ds05$V12, sep=""))
ds09$V13<-as.factor(paste(ds09$V11,ds09$V12, sep=""))


res_err00<-do.call(rbind,by(ds00, ds00$V1, function(x) expected_error_rate(x)))
res_err01<-do.call(rbind,by(ds01, ds01$V1, function(x) expected_error_rate(x)))
res_err05<-do.call(rbind,by(ds05, ds05$V1, function(x) expected_error_rate(x)))
res_err09<-do.call(rbind,by(ds09, ds09$V1, function(x) expected_error_rate(x)))

pdf(file = "Dstat_error",width = 8,height = 8)
par(mfrow=c(2,2))
boxplot(res_err00, main = 'Dstat error (extinction rate = 0.0)')
boxplot(res_err01, main = 'Dstat error (extinction rate = 0.1)')
boxplot(res_err05, main = 'Dstat error (extinction rate = 0.5)')
boxplot(res_err09, main = 'Dstat error (extinction rate = 0.9)')
dev.off()

d300=read.table('../D3/D3_00',h=F)
d301=read.table('../D3/D3_01',h=F)
d305=read.table('../D3/D3_05',h=F)
d309=read.table('../D3/D3_09',h=F)
d300$V8<-as.factor(paste(d300$V6,d300$V7, sep=""))
d301$V8<-as.factor(paste(d301$V6,d301$V7, sep=""))
d305$V8<-as.factor(paste(d305$V6,d305$V7, sep=""))
d309$V8<-as.factor(paste(d309$V6,d309$V7, sep=""))



pvalue<-function(lD3,D3) {
  pvalue<-length(lD3[lD3>D3])/length(lD3)
  return(pvalue)
}

expected3_error_rate<-function(temp){
  n_H0<-nrow(temp[temp$V8 == "P3P2" |
  temp$V8 == "P2P3"|
  temp$V8 == "P3P1"|
  temp$V8 == "P1P3",])
  n_alt<-nrow(temp[temp$V8 == "OP1"|temp$V8 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)

  nulD3<-temp[temp$V8 == "P2P2"|
  temp$V8 == "P1P1"|
  temp$V8 == "P3P3"|
  temp$V8 == "N1N1"|
  temp$V8 == "N1P3"|
  temp$V8 == "N1O"|
  temp$V8 == "P1O"|
  temp$V8 == "P1N1"|
  temp$V8 == "P2N1"|
  temp$V8 == "P3N1"|
  temp$V8 == "P2O"|
  temp$V8 == "P3O"|
  temp$V8 == "OO"|
  temp$V8 == "OP3",]

  temp$V9<-unlist(lapply(temp$V5,function(x) pvalue(abs(nulD3$V5),abs(x))))
  temp$V10<-(p.adjust(temp$V9,"fdr"))
  signi_D<-temp[temp$V10<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V8 == "P3P2" |
  signi_D$V8 == "P2P3"|
  signi_D$V8 == "P3P1"|
  signi_D$V8 == "P1P3",])
  err_all<- 1 - (true_pos/nrow(signi_D))

  alt=nrow(signi_D[signi_D$V8 == "OP2" |
  signi_D$V8 == "OP1"|
  signi_D$V8 == "N2P1"|
  signi_D$V8 == "N2P2",])
  err_alt=alt/nrow(signi_D)
  data=c("expected_error" = exp_err,"observed_error" = err_alt)
  return(data)
}

temp=d300[d300$V1=='test2000',]

res3_err00<-do.call(rbind,by(d300, d300$V1, function(x) expected3_error_rate(x)))
res3_err01<-do.call(rbind,by(d301, d301$V1, function(x) expected3_error_rate(x)))
res3_err05<-do.call(rbind,by(d305, d305$V1, function(x) expected3_error_rate(x)))
res3_err09<-do.call(rbind,by(d309, d309$V1, function(x) expected3_error_rate(x)))

pdf(file = "D3_error",width = 8,height = 8)
par(mfrow=c(2,2))
boxplot(res3_err00, main = 'D3 error (extinction rate = 0.0)')
boxplot(res3_err01, main = 'D3 error (extinction rate = 0.1)')
boxplot(res3_err05, main = 'D3 error (extinction rate = 0.5)')
boxplot(res3_err09, main = 'D3 error (extinction rate = 0.9)')
dev.off()


test<-"test2000"
temp<-ds00[ds00$V1==test,]
temp<-temp[order(temp$V13),]
table(temp$V13)
plot(temp$V8, pch=19,cex=0.5, col = "red")
abline(h=median(temp$V8), col='blue')
abline(h=0,col="red")
