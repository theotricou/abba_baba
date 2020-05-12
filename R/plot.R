#!/usr/bin/env Rscript
# Theo

#shell

for i in `seq 2678 1 2678`; do
  sed "s/aaaa/$i/g" run_slurm.sh > temp
  # sbatch temp
  # rm temp
done

rm -rf D*
for j in sim_[0159];do
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

################################################################################
################################################################################
################################################################################


# Dz error rate
error_Ds<-function(temp){
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|temp$V14 == "P1P3",])
  n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)
  temp$V15<-(p.adjust(temp$V9,"fdr"))
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V15<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  alt=nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  err_alt=alt/(true_pos+alt)
  sen_H0<-true_pos/n_H0
  sen_Ha<-alt/n_alt
  err_all=true_pos/nrow(signi_D)
  data=c("s" = levels(dataD$V1)[temp$V1[1]],
  # "t" = levels(dataD[dataD$V1 == levels(dataD$V1)[temp$V1[1]], 2])[temp$V2[1]],
  "n" =n_H0+n_alt , "err_H0"=sen_H0,"err_Ha"=sen_Ha,"err_D"=err_all,"err_exp" = exp_err,"err_obs" = err_alt)
  return(data)
}

pvalue<-function(lD3,D3) {
  pvalue<-length(lD3[lD3>D3])/length(lD3)
  return(pvalue)
}

error_D3<-function(temp){
  n_H0<-nrow(temp[temp$V9 == "P3P2"|temp$V9 == "P2P3"|temp$V9 == "P3P1"|temp$V9 == "P1P3",])
  n_alt<-nrow(temp[temp$V9 == "OP1"|temp$V9 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)
  nulD3<-temp[temp$V9 == "P2P2"|temp$V9 == "P1P1"|temp$V9 == "P3P3"|temp$V9 == "N1N1"|
    temp$V9 == "N1P3"|temp$V9 == "N1O"|temp$V9 == "P1O"|temp$V9 == "P1N1"|
    temp$V9 == "P2N1"|temp$V9 == "P3N1"|temp$V9 == "P2O"|temp$V9 == "P3O"|
    temp$V9 == "OO"|temp$V9 == "OP3",]
  temp$V10<-unlist(lapply(temp$V6,function(x) pvalue(abs(nulD3$V6),abs(x))))
  temp$V11<-p.adjust(temp$V10,"fdr")
  # temp$V11<-temp$V10
  signi_D<-temp[temp$V10<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V9 == "P3P2"|signi_D$V9 == "P2P3"|
    signi_D$V9 == "P3P1"|signi_D$V9 == "P1P3",])
  alt=nrow(signi_D[signi_D$V9 == "OP2"|signi_D$V9 == "OP1",])
  err_alt=alt/(true_pos+alt)
  sen_H0<-true_pos/n_H0
  sen_Ha<-alt/n_alt
  err_all=true_pos/nrow(signi_D)
  data=c("s" = levels(dataD3$V1)[temp$V1[1]],"n"=n_H0+n_alt, "err_H0"=sen_H0,"err_Ha"=sen_Ha,"err_D"=err_all,"err_exp" = exp_err,"err_obs" = err_alt)
  return(data)
}

dataD<-read.table('Dstat',h=F)
dataD$V14<-as.factor(paste(dataD$V12,dataD$V13, sep=""))
resD<-as.data.frame(do.call(rbind,by(dataD, dataD[,1:2], function(x) error_Ds(x))))
library(reshape2)
resD$s<-as.factor(resD$s)
temp<-resD[,c(1,6,7)]
colnames(temp)<-c('ex','obs','exp')
temp$obs<-as.numeric(as.character(temp$obs))
temp$exp<-as.numeric(as.character(temp$exp))
temp = (melt(temp, id=("ex")))


temp<-resD[,c(5,6)]
colnames(temp)<-c('obs','exp')
temp$obs<-as.numeric(as.character(temp$obs))
temp$exp<-as.numeric(as.character(temp$exp))
temp = (melt(temp, id=("ex")))


pdf(file = "Dstat_error.pdf", width = 8, height = 8)
boxplot(temp$value~temp$variable+temp$ex,las = 1,
  main="D stat error (expected vs observed) by extinction rate",
  col = rep(c("white", "grey"), each = 1),
    names=c("exp\n0.0","obs\n0.0","exp\n0.1","obs\n0.1","exp\n0.5","obs\n0.5",
    "exp\n0.9","obs\n0.9"),cex.axis = 0.75,
   yaxt="n" )
axis(side = 2, at = pretty(0:1,10))
dev.off()



dataD3<-read.table('D3',h=F)
dataD3$V9<-as.factor(paste(dataD3$V7,dataD3$V8, sep=""))
res3<-as.data.frame(do.call(rbind,by(dataD3, dataD3[,1:2], function(x) error_D3(x))))

library(reshape2)
res3$s<-as.factor(resD$s)
temp<-res3[,c(1,6,7)]
colnames(temp)<-c('ex','obs','exp')
temp$obs<-as.numeric(as.character(temp$obs))
temp$exp<-as.numeric(as.character(temp$exp))
temp = (melt(temp, id=("ex")))


pdf(file = "D3_error.pdf", width = 8, height = 8)
boxplot(temp$value~temp$variable+temp$ex,las = 1,
  main="D3 error (expected vs observed) by extinction rate",
  col = rep(c("white", "grey"), each = 1),
    names=c("exp\n0.0","obs\n0.0","exp\n0.1","obs\n0.1","exp\n0.5","obs\n0.5",
    "exp\n0.9","obs\n0.9"),cex.axis = 0.75,
   yaxt="n" )
axis(side = 2, at = pretty(0:1,10))
dev.off()

# Ds error by distance P3P4

Ds_error_by_dist<-function(temp){
  ### expected
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|temp$V14 == "P1P3",])
  n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)

  ### observed
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V13<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  alt<-nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  err_alt<-alt/(true_pos+alt)
  return(exp_err)
}

truc<-function(t) {
  t$V13<-p.adjust(t$V9,"fdr")
  res<-as.data.frame(rbind(by(t, t$V15, function(x) Ds_error_by_dist(x))))
  return(res)
}


par(mfrow=c(2,2))


block=100000
pdf(file = "Dstat_dist_err.pdf", width = 8, height = 8)
par(mfrow=c(2,2))
dataDist<-read.table('sim_0/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
boxplot(ed, main= "ext = 0.0")

dataDist<-read.table('sim_1/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
boxplot(ed, main= "ext = 0.1")

dataDist<-read.table('sim_5/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
boxplot(ed, main= "ext = 0.5")

dataDist<-read.table('sim_9/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
boxplot(ed, main= "ext = 0.9")

dev.off()




# test graphic

dataDist<-read.table('sim_0/Dstat',h=F)
dataDist$V13<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
tempDs<-dataDist[order(dataDist$V13,dataDist$V8),]
plot(tempDs$V8)
table(tempDs$V13)



dataD3<-read.table('sim_0/D3',h=F)
dataD3$V8<-as.factor(paste(dataD3$V6,dataD3$V7, sep=""))
tempD3<-dataD3[order(dataD3$V8,dataD3$V5),]
plot(tempD3$V5)
table(tempD3$V8)




mean_ev_D3<-as.data.frame(cbind(by(tempD3, tempD3$V8, function(x) mean(x[,5]))))
mean_ev_D3$V2<-rownames(mean_ev_D3)

pdf(file = "D3_mean_ev.pdf", width = 16, height = 8)
plot(mean_ev_D3$V1,xaxt="n")
axis(side = 1, at = seq(1:25),labels = mean_ev_D3$V2, cex.axis = 0.75)
dev.off()


mean_ev_Ds<-as.data.frame(cbind(by(dataDist, dataDist$V13, function(x) mean(x[,8]))))
mean_ev_Ds$V2<-rownames(mean_ev_Ds)

pdf(file = "Ds_mean_ev.pdf", width = 16, height = 8)
plot(mean_ev_Ds$V1,xaxt="n")
axis(side = 1, at = seq(1:nrow(mean_ev_Ds)),labels = mean_ev_Ds$V2, cex.axis = 0.45)
dev.off()

pdf(file = "Ds_mean_ev_V2.pdf", width = 16, height = 8)
boxplot(V8~V13, data=dataDist,cex.axis = 0.45,outline=FALSE)
dev.off()

pdf(file = "D3_mean_ev_V2.pdf", width = 16, height = 8)
boxplot(V5~V8, data=dataD3,cex.axis = 0.75, outline=FALSE)
dev.off()
