# Dz error rate
error_Ds<-function(temp){
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|temp$V14 == "P1P3",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
  exp_err<-n_alt/(n_alt+n_H0)
  temp$V15<-(p.adjust(temp$V9,"fdr"))
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V15<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  # alt=nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  alt=nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
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
mean_err<-rbind(by(temp, temp[,c(2,1)], function(x) mean(na.omit(x$value))))




resD2<-resD[resD$err_exp != 0,]
temp<-resD2[,c(1,6,7)]
colnames(temp)<-c('ex','obs','exp')
temp$obs<-as.numeric(as.character(temp$obs))
temp$exp<-as.numeric(as.character(temp$exp))
temp = (melt(temp, id=("ex")))
mean_err<-rbind(by(temp, temp[,c(2,1)], function(x) mean(na.omit(x$value))))




pdf(file = "Dstat_error_no_out.pdf", width = 8, height = 8)
boxplot(temp$value~temp$variable+temp$ex,las = 1,
  main="D stat error (expected vs observed) by extinction rate",
  col = rep(c("white", "grey"), each = 1),
    names=c("exp\n0.0","obs\n0.0","exp\n0.1","obs\n0.1","exp\n0.5","obs\n0.5",
    "exp\n0.9","obs\n0.9"),cex.axis = 0.75,
   yaxt="n" )
axis(side = 2, at = pretty(0:1,10))
points(mean_err[c(2,1,4,3,6,5,8,7)],col="red",pch=18, cex = 2)
dev.off()




pdf(file = "one_tree_error_v2.pdf", width = 8, height = 8)
boxplot(temp$value~temp$variable+temp$ex,las = 1,
  main="D stat error (expected vs observed) by extinction rate",
  col = rep(c("white", "grey"), each = 1),
  cex.axis = 0.75, yaxt="n" )
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
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
  n_alt<-nrow(temp[temp$V14 == "OP1"|temp$V14 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)

  ### observed
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V13<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  # alt<-nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  # alt<-nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  alt<-nrow(signi_D[signi_D$V14 == "OP1"|signi_D$V14 == "OP2",])
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
block=200000
dataDist<-read.table('sim_0/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.0")
points(m_ed, col='red')
dataDist<-read.table('sim_1/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.1")
points(m_ed, col='red')
dataDist<-read.table('sim_5/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.5")
points(m_ed, col='red')
dataDist<-read.table('sim_9/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
dataDist$V15 = cut(dataDist$V10,seq(0,1000000,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.9")
points(m_ed, col='red')

dev.off()



# Ds error by RELATIVE distance P3P4

Ds_error_by_dist<-function(temp){
  ### expected
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|
    temp$V14 == "P1P3",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
  n_alt<-nrow(temp[temp$V14 == "OP1"|temp$V14 == "OP2",])
  exp_err<-n_alt/(n_alt+n_H0)

  ### observed
  # temp$V14<-temp$V9
  signi_D<-temp[temp$V13<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|
    signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  # alt<-nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  # alt<-nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  alt<-nrow(signi_D[signi_D$V14 == "OP1"|signi_D$V14 == "OP2",])
  err_alt<-alt/(true_pos+alt)
  return(exp_err)
}

truc<-function(t) {
  t$V13<-p.adjust(t$V9,"fdr")
  res<-as.data.frame(rbind(by(t, t$V15, function(x) Ds_error_by_dist(x))))
  return(res)
}

relative<-function(x){
  res<-x$V10/max(x$V10)
  return(res)
}

block=0.1
dataDist<-read.table('sim_5/Dstat',h=F)
dataDist$V14<-as.factor(paste(dataDist$V11,dataDist$V12, sep=""))
# RD calc
dataDist$RD<-unlist(by(dataDist, dataDist$V1, function(x) relative(x)))
dataDist$V15 = cut(dataDist$RD,seq(0,1,block))
dataDist$V15<-as.factor(dataDist$V15)
ed<-do.call(rbind,by(dataDist, dataDist$V1, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.0")
points(m_ed, col='red')



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









library('apTreeshape')


t1=read.tree('spe09')
extant<-grep("l>",t1$tip.label)
t1<-as.treeshape(keep.tip(t1, t1$tip.label[extant]))
colless(t1)


colless(t1)











############
#Sampling effect
#One tree, x transfers events, 10 random sampling rate



set.seed(123456789)
data<-read.table('oneDstat2',h=F)
data$V14<-as.factor(paste(data$V12,data$V13, sep=""))
ind<-unique(c(unique(data$V3),unique(data$V4),unique(data$V5),unique(data$V6)))





lvl<-seq(4,length(ind),((length(ind)-4)/10))

lapply(lvl, function(x) sampling_error(x))

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

sampling_error<-function(d){
  temp<-sample(ind ,d ,replace = FALSE)
  sub<-subset(d, V3 %in% temp & V4 %in% temp & V5 %in% temp & V6 %in% temp))
  resD<-as.data.frame(do.call(rbind,by(dataD, dataD[,1:2], function(x) error_Ds(x))))
}





data<-read.table('tempDstat',h=F)
data$V14<-as.factor(paste(data$V12,data$V13, sep=""))
ind<-unique(c(unique(data$V3),unique(data$V4),unique(data$V5),unique(data$V6)))
data$V15<-(p.adjust(data$V9,"fdr"))
lvl<-seq(4,length(ind),4)

samp<-function(x){
ll=c()
for (i in 1:200){
  # set.seed(i)
  temp<-sample(ind ,x,replace = FALSE)
  sub<-subset(data, V3 %in% temp & V4 %in% temp & V5 %in% temp & V6 %in% temp)
  signi_D<-sub[sub$V15<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|
    signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  alt=nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|
    signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  err_alt=alt/(true_pos+alt)
  ll[i]<-err_alt
}
return(ll)
}

bb=do.call(rbind, lapply(lvl, function(x) samp(x)))
boxplot(t(bb), ylim = c(0,1), xaxt="n")
axis(side = 1, at = seq(1:length(lvl)),labels = lvl, cex.axis = 1)
means <- apply(bb,1, function(x) mean(na.omit(x)))
points(means,col="red",pch=18, cex = 2)




pdf(file = "sampling_00.pdf", width = 8, height = 8)
boxplot(t(bb), ylim = c(0,1), xaxt="n")
axis(side = 1, at = seq(1:length(lvl)),labels = lvl, cex.axis = 1)
means <- apply(bb,1, function(x) mean(na.omit(x)))
points(means,col="red",pch=18, cex = 2)
dev.off()











for i in test*; do
  sed "1d" $i/data.txt | sed "s/^/$i\t/" >> oneDstat
  sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> oneD3
done

sed "s/^/09\t/" oneDstat > oneDstat2




error_Ds<-function(temp){
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|temp$V14 == "P1P3",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
  exp_err<-n_alt/(n_alt+n_H0)
  temp$V15<-(p.adjust(temp$V9,"fdr"))
  # temp$V15<-temp$V9
  signi_D<-temp[temp$V15<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  # alt=nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  alt=nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  err_alt=alt/(true_pos+alt)
  sen_H0<-true_pos/n_H0
  sen_Ha<-alt/n_alt
  err_all=true_pos/nrow(signi_D)
  data=c("s" = levels(dataD$V1)[temp$V1[1]],
  # "t" = levels(dataD[dataD$V1 == levels(dataD$V1)[temp$V1[1]], 2])[temp$V2[1]],
  "n" =n_H0+n_alt , "err_H0"=sen_H0,"err_Ha"=sen_Ha,"err_D"=err_all,"err_exp" = exp_err,"err_obs" = err_alt)
  return(data)
}




dataD<-read.table('oneDstat2',h=F)
dataD$V14<-as.factor(paste(dataD$V12,dataD$V13, sep=""))
resD<-as.data.frame(do.call(rbind,by(dataD, dataD[,1:2], function(x) error_Ds(x))))

par(mfrow=c(1,1))

pdf(file = "error_1000_tr.pdf", width = 8, height = 8)
boxplot(resD[5:6], ylim = c(0,1))
means <- apply(resD[5:6],2, function(x) mean(na.omit(x)))
points(means,col="red",pch=18, cex = 2)
dev.off()


err_by_quat<-function(temp){
  n_H0<-nrow(temp[temp$V14 == "P3P2"|temp$V14 == "P2P3"|temp$V14 == "P3P1"|temp$V14 == "P1P3",])
  # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
  n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
  exp_err<-n_alt/(n_alt+n_H0)
  temp$V15<-(p.adjust(temp$V9,"fdr"))
  # temp$V15<-temp$V9
  signi_D<-temp[temp$V15<=0.05,]
  true_pos<-nrow(signi_D[signi_D$V14 == "P3P2"|signi_D$V14 == "P2P3"|signi_D$V14 == "P3P1"|signi_D$V14 == "P1P3",])
  # alt=nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  alt=nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
  err_alt=alt/(true_pos+alt)
  data=c("err_exp" = exp_err,"err_obs" = err_alt)
}


Q<-dataD[dataD$V2=="test2501",c(3,4,5,6)]
# x<-Q[1671,]
res_err_quat<-t(as.data.frame(apply(Q,1,function(x) c(x, err_by_quat(subset(dataD, V3 %in% x[1] & V4 %in% x[2] & V5 %in% x[3] & V6 %in% x[4]))))))


pdf(file = "quat_err.pdf", width = 16, height = 8)
par(mfrow=c(1,2))
plot(temp[,c(1,2)])
lines(x=1,y=1, col = 'red')
plot(temp[,c(1,3)])
dev.off()


pdf(file = "proba.pdf", width = 8, height = 8)
boxplot(result, ylim = c(0,1))
dev.off()
