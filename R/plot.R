#!/usr/bin/env Rscript
# Theo

library('ape')
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

expected_non_nul3<-which(!d3$Recip == "O" &
!d3$Recip == "N1"&
!(d3$Donor == "P1" & d3$Recip== "P2")&
!(d3$Donor == "P2" & d3$Recip== "P1"))

d3_non_nul<-d3[expected_non_nul3,]
d3_nul<-d3[setdiff(seq(1,nrow(d3)), expected_non_nul3),]

ks.test(d3_non_nul[1,]$D3, d3_nul$D3)


library("fdrtool")
ss = seq(1,length(d3$D3))
W1 <- lm(d3$D3 ~ ss)
CR <- (d3$D3-W1$fitted.values)/summary(W1)$sigma
# mydata$V6 <- (mydata$V3-W1$fitted.values)/summary(W1)$sigma
normCR <- 1-pnorm(CR)
adjustnormCR <- p.adjust(normCR, "fdr")
d3[as.numeric(as.character(names(adjustnormCR[adjustnormCR <= 0.05]))),]

plot(d3$D3)

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
