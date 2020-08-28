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

# Ds_error_by_dist<-function(temp){
#   ### expected
#   n_H0<-nrow(temp[temp$EV == "P3P2"|temp$EV == "P2P3"|temp$EV == "P3P1"|
#     temp$EV == "P1P3",])
#   # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "OP1"|temp$V14 == "N2P1"|temp$V14 == "OP2",])
#   # n_alt<-nrow(temp[temp$V14 == "N2P2"|temp$V14 == "N2P1",])
#   n_alt<-nrow(temp[temp$EV == "N2P1"|temp$EV == "N2P2",])
#   exp_err<-n_alt/(n_alt+n_H0)
#
#   ### observed
#   # temp$V14<-temp$V9
#   signi_D<-temp[temp$Adjust<=0.05,]
#   true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
#     signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
#   # alt<-nrow(signi_D[signi_D$V14 == "OP2"|signi_D$V14 == "OP1"|signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
#   # alt<-nrow(signi_D[signi_D$V14 == "N2P1"|signi_D$V14 == "N2P2",])
#   alt<-nrow(signi_D[signi_D$EV == "N2P1"|signi_D$EV == "N2P2",])
#   obs_err<-alt/(true_pos+alt)
#   return(obs_err)
# }

adjus<-function(t) {
  return(p.adjust(t$Pvalue,"fdr"))
}

truc<-function(t) {
  res<-as.data.frame(rbind(by(t, t$Block, function(x) Ds_error_by_dist(x))))
  return(res)
}

dist_relative<-function(x){
  N2<-x$dP1P4 - x$dP1P3
  res<-N2/max(N2)
  return(res)
}

Ds_error_by_dist<-function(temp){
  ### expected
  n_H0<-nrow(temp[temp$EV == "P3P2"|temp$EV == "P2P3"|temp$EV == "P3P1"|
    temp$EV == "P1P3",])
  n_alt<-nrow(temp[temp$EV == "N2P1"|temp$EV == "N2P2",])
  exp_err<-n_alt/(n_alt+n_H0)

  ### observed
  signi_D<-temp[temp$Adjust<=0.05,]
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
  alt<-nrow(signi_D[signi_D$EV == "N2P1"|signi_D$EV == "N2P2",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}

out<-function(temp){
  ### observed
  signi_D<-temp[temp$Adjust<=0.05,]
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
  alt<-nrow(signi_D[signi_D$EV == "OP1"|signi_D$EV == "OP2",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}

truc2<-function(t) {
  res<-as.data.frame(rbind(by(t, t$Block, function(x) out(x))))
  return(res)
}


# par(mfrow=c(1,3))


dataDist<-read.table('new_05/data_D.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
# RD calc
dataDist$RD<-unlist(by(dataDist, dataDist$test, function(x) dist_relative(x)))
block=0.1
dataDist$Block = as.factor(cut(dataDist$RD,seq(0,1,block)))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))
ed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))
oed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc2(x)))
om_ed<-apply(oed,2,function(x) mean(na.omit(x)))
boxplot(ed, main= "ext = 0.5")
points(m_ed, col='red')
points(om_ed, col='blue')



block=100000
dataDist$Block = as.factor(cut(dataDist$dist,seq(0,1000000,block)))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))
ed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))










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











################################################################################


library("ape")
library("phangorn")
library("parallel")

colless<-function(tr) {
    leftright4all<-t(sapply((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)),function(x,tree) Children(tree,x), tree=tr))
    leftright4all.size<-apply(leftright4all,2,function(x,tree) unlist(lapply(Descendants(tree, x, "tips"), length)),tree=tr)
    index<-sum(abs(apply(leftright4all.size,1,diff)))
    return(index)
}

collesser<-function(tr){
  extantsampled<-tr$tip.label[grep("l>",tr$tip.label)]
  extant.tr<-keep.tip(tr, extantsampled)
  res<-c("c_col"=as.integer(colless(tr)),"e_col"=as.integer(colless(extant.tr)))
  return(res)
}

err_D<-function(temp){
  D_H0<-length(which(temp$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  D_H1<-length(which(temp$EV %in% c('N2P1','N2P2')))
  D_H2<-length(which(temp$EV %in% c('OP1','OP2')))
  D_H3<-length(which(temp$EV %in% c('P4P1','P4P2')))
  D_err_exp_H1<-D_H1/(D_H0+D_H1)
  D_err_exp_H2<-D_H2/(D_H0+D_H2)
  D_err_exp_H3<-D_H3/(D_H0+D_H3)
  D_err_exp_H123<-(D_H1+D_H2+D_H3)/(D_H0+D_H1+D_H2+D_H3)
  stemp<-temp[temp$FDR<=0.05,]
  sD_H0<-length(which(stemp$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  sD_H1<-length(which(stemp$EV %in% c('N2P1','N2P2')))
  sD_H2<-length(which(stemp$EV %in% c('OP1','OP2')))
  sD_H3<-length(which(stemp$EV %in% c('P4P1','P4P2')))
  sD_err_exp_H1<-sD_H1/(sD_H0+sD_H1)
  sD_err_exp_H2<-sD_H2/(sD_H0+sD_H2)
  sD_err_exp_H3<-sD_H3/(sD_H0+sD_H3)
  sD_err_exp_H123<-(sD_H1+sD_H2+sD_H3)/(sD_H0+sD_H1+sD_H2+sD_H3)
  res<-c("Derr_exp"=D_err_exp_H1,"DH2err_exp"=D_err_exp_H2,"DH3err_exp"=D_err_exp_H3,"DH123err_exp"=D_err_exp_H123,
    "Derr_obs"=sD_err_exp_H1,"DH2err_obs"=sD_err_exp_H2,"DH3err_obs"=sD_err_exp_H3,"DH123err_obs"=sD_err_exp_H123)
  return(res)
}

err_D3<-function(temp){
  sub<-temp[temp$EV %in% c("P1P1","P2P2","P3P3","N1N1","OO","OP3","P1O","P2O","P3O","P1N1","P2N1","P3N1"),]
  basic_pvalue<-function(D3) {
    pvalue<-length(which((abs(sub$D3) + (7*sd(abs(sub$D3)))) >= abs(D3)))/nrow(sub)
    return(pvalue)
  }
  temp$Pvalue<-unlist(lapply(temp$D3,function(x) basic_pvalue(x)))
  D3_H0<-length(which(temp$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  D3_H1<-length(which(temp$EV %in% c('OP1','OP2')))
  D3_err_exp_H1<-D3_H1/(D3_H0+D3_H1)
  stemp<-temp[which(p.adjust(temp$Pvalue,"fdr") <= 0.05),]
  stemp<-temp[temp$Pvalue <= 0.05,]
  sD3_H0<-length(which(stemp$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  sD3_H1<-length(which(stemp$EV %in% c('OP1','OP2')))
  sD3_err_exp_H1<-sD3_H1/(sD3_H0+sD3_H1)
  res<-c("D3err_exp"=D3_err_exp_H1,"D3err_obs"=sD3_err_exp_H1)
  return(res)
}

validateandreorder<-function(arr, dist) {
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  }
}
getquatuors<-function(tr) {
  dist<-cophenetic(tr)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}

overlap<-function(ab,cd) {
    return(!(ab[2]<cd[1] | cd[2]<ab[1]))
}
getstartstopofedges<-function(mat,xx) {
    if (!is.null(mat)) {
        return(t(apply(mat,1,function(x,xx) c(x,xx[x]),xx=xx)))
    }
    else redataDist<-read.table('new_05/data_D.txt',h=T)
turn(NULL)
}

plottreeandinfo<-function(tr, quat, returnxy=TRUE, plot=TRUE) {
    plot(tr, plot=plot)
    xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (returnxy) return(xy)
}

probas.donor.receiver<-function(quat, tr, xy) { #xy are x and y coordinates in the plot
    p1<-quat[1]
    p2<-quat[2]
    p3<-quat[3]
    p4<-quat[4]
    p1p2<-mrca.phylo(tr, c(p1,p2))
    p1nb<-which(tr$tip.label==p1)
    p2nb<-which(tr$tip.label==p2)
    desc<-sapply(Children(tr, p1p2), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p2nb, desc[[1]])) desc<-rev(desc) ##so that p1 remains in first position.
    bl1.donor<-sum(tr$edge.length[match(desc[[1]],tr$edge[,2])])
    bl2.donor<-sum(tr$edge.length[match(desc[[2]],tr$edge[,2])])
    P1DONOR<-tr$edge[match(desc[[1]],tr$edge[,2]),, drop=FALSE]
    P2DONOR<-tr$edge[match(desc[[2]],tr$edge[,2]),,drop=FALSE]
    p1.up<-c(p1nb, Ancestors(tr, p1nb))
    p1.up.small<-p1.up[1:(which(p1.up==p1p2)-1)]
    p2.up<-c(p2nb, Ancestors(tr, which(tr$tip.label==p2)))
    p2.up.small<-p2.up[1:(which(p2.up==p1p2)-1)]
    bl1.receiver = sum(tr$edge.length[match(p1.up.small,tr$edge[,2])])
    bl2.receiver = sum(tr$edge.length[match(p2.up.small,tr$edge[,2])])
    P1RECEIVER<-tr$edge[match(p1.up.small,tr$edge[,2]),, drop=FALSE]
    P2RECEIVER<-tr$edge[match(p2.up.small,tr$edge[,2]),, drop=FALSE]
    p1p3<-mrca.phylo(tr, c(p1,p3)) #same than p2p3
    p3nb<-which(tr$tip.label==p3)
    desc1.3<-sapply(Children(tr, p1p3), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p1nb, desc1.3[[1]])) desc1.3<-rev(desc1.3)
    desc1.3<-desc1.3[[1]]
    subedge<-tr$edge[match(desc1.3,tr$edge[,2]),, drop=FALSE]
    xp1p2<-xy$xx[p1p2] #time of p1p2 speciation. p3 donor or receiver must be more recent.
    xp1<-xy$xx[p1nb]
    subedgetime<-cbind(subedge, t(apply(subedge,1,function(x,xx) xx[x], xx=xy$xx)))
    subedgetime<-subedgetime[(subedgetime[,3]>xp1p2)|(subedgetime[,4]>xp1p2),, drop=FALSE] ##select those that have one side at least in the range encounterd by p1p2
    subedgetime[subedgetime[,3]<xp1p2,3]<-xp1p2 #change start time of the edges traversing the line to only count possible bl.
    bl3.donor<-sum(subedgetime[,4]-subedgetime[,3])
    P3DONOR<-subedge
    p3.up<-c(p3nb, Ancestors(tr, p3nb))
    p3.up.subedgetime<-match(p3.up, subedgetime[,2])
    subedgetime.up<-subedgetime[p3.up.subedgetime[!is.na(p3.up.subedgetime)],, drop=FALSE]
    bl3.receiver<-sum(subedgetime.up[,4]-subedgetime.up[,3])
    P3RECEIVER<-subedgetime.up[,1:2, drop=FALSE]
    p1p4<-mrca.phylo(tr, c(p1,p4))
    p4nb<-which(tr$tip.label==p4)
    possible.p4.nodes<-intersect(Descendants(tr, p1p4, "all"), Ancestors(tr, p1p3))
    if (length(possible.p4.nodes)==0) {
        bl4.donor<-0
        P4DONOR<-NULL
    }
    else {
        p4.to.check<-setdiff(unique(unlist(Descendants(tr, possible.p4.nodes,"all"))), Descendants(tr, p1p3,"all"))
        subedge.p4<-tr$edge[match(p4.to.check,tr$edge[,2]),]
        subedgetime.p4<-cbind(subedge.p4, t(apply(subedge.p4,1,function(x,xx) xx[x], xx=xy$xx)))
        subedgetime.p4<-subedgetime.p4[(subedgetime.p4[,3]>xp1p2)|(subedgetime.p4[,4]>xp1p2),,drop=FALSE] ##select those that have one side at least in the range encounterd by p1p2
        subedgetime.p4[subedgetime.p4[,3]<xp1p2,3]<-xp1p2 #change start time of the edges traversing the line to only count possible bl.
        bl4.donor<-sum(subedgetime.p4[,4]-subedgetime.p4[,3])
        P4DONOR<-subedgetime.p4[,1:2,drop=FALSE]
        if (nrow(P4DONOR)==0) P4DONOR<-NULL
    }
    relativedist2outgroup<-(xy$xx[p1p3]-xy$xx[p1p4])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top3<-(xy$xx[p1nb]-xy$xx[p1p3])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top2<-(xy$xx[p1nb]-xy$xx[p1p2])/(xy$xx[p1nb]-xy$xx[p1p4])
    treesize<-sum(tr$edge.length)
    pP1d<-bl1.donor/treesize
    pP2d<-bl2.donor/treesize
    pP3d<-bl3.donor/treesize
    pP4d<-bl4.donor/treesize
    pP1r<-bl1.receiver/treesize
    pP2r<-bl2.receiver/treesize
    pP3r<-bl3.receiver/treesize
    DONOR.RECEIVER<-list(P1DONOR=P1DONOR, P2DONOR=P2DONOR,P3DONOR=P3DONOR,P4DONOR=P4DONOR,P1RECEIVER=P1RECEIVER,P2RECEIVER=P2RECEIVER,P3RECEIVER=P3RECEIVER)
    DONOR.RECEIVER<-lapply(DONOR.RECEIVER, getstartstopofedges, xx=xy$xx)
    getpossibletransfers<-function(mat1,mat2) {
        if (!is.null(mat1)&!is.null(mat2)) {
            TFmat<-apply(mat1[,3:4, drop=FALSE],1,function(x,m2) apply(m2[,3:4, drop=FALSE], 1, overlap,cd=x), m2=mat2)
            if (is.null(dim(TFmat))) TFmat<-t(t(TFmat))
            else TFmat<-t(TFmat)
            ##in the matrix, rows are rows of mat1 and cols are cols of
            matchingedges<-which(TFmat, arr.ind=T)
            ##and for each match, we recover its condensed name.
            res<-apply(matchingedges,1,function(x,m1,m2) paste(paste(m1[x[1],1:2],collapse="-"),paste(m2[x[2],1:2],collapse="-"),sep="->"), m1=mat1, m2=mat2)
            return(res)
        }
    }
    P1P3<-getpossibletransfers(DONOR.RECEIVER$P1DONOR, DONOR.RECEIVER$P3RECEIVER)
    P2P3<-getpossibletransfers(DONOR.RECEIVER$P2DONOR, DONOR.RECEIVER$P3RECEIVER)
    P3P1<-getpossibletransfers(DONOR.RECEIVER$P3DONOR, DONOR.RECEIVER$P1RECEIVER)
    P3P2<-getpossibletransfers(DONOR.RECEIVER$P3DONOR, DONOR.RECEIVER$P2RECEIVER)
    P4P1<-getpossibletransfers(DONOR.RECEIVER$P4DONOR, DONOR.RECEIVER$P1RECEIVER)
    P4P2<-getpossibletransfers(DONOR.RECEIVER$P4DONOR, DONOR.RECEIVER$P2RECEIVER)
    TRANSFERSINVOLVED<-c(P1P3,P2P3,P3P1,P3P2,P4P1,P4P2)
    return(list(pP1d = pP1d,pP2d = pP2d,pP3d = pP3d,pP4d = pP4d,pP1r = pP1r,pP2r = pP2r,pP3r = pP3r, relativedist2outgroup=relativedist2outgroup, relativedistp1top3=relativedistp1top3,relativedistp1top2=relativedistp1top2, transfersinvolved=TRANSFERSINVOLVED))
}

ProbaFalsePositive<-function(P) {
    proba.true.positive<-(P$pP1d*P$pP3r)+(P$pP3d*P$pP1r)+(P$pP2d*P$pP3r)+(P$pP3d*P$pP2r)
    proba.false.positive<-(P$pP4d*P$pP2r)+(P$pP4d*P$pP1r)
    return(list(tp=proba.true.positive, fp=proba.false.positive))
}

tree_proba_err<-function(tr){
  extantsampled<-tr$tip.label[grep("l>",tr$tip.label)]
  extant.tr<-keep.tip(tr, extantsampled)
  quatuors <- getquatuors(extant.tr)
  xy<-plottreeandinfo(tr, quatuors[1,], returnxy=TRUE, plot=T)
  PROBAS<-apply(quatuors, 1, probas.donor.receiver, tr=tr, xy=xy)
  PPP<-lapply(PROBAS, ProbaFalsePositive)
  RES<-do.call(rbind, lapply(PPP, unlist))
  result<-RES[,2]/(RES[,1]+RES[,2])
  data<-c("proba"=mean(na.omit(result)))
  return(data)
}

block_err<-function(temp){
  block=1000000/6
  lvl<-levels(cut(1000000, seq(0,1000000,block)))
  temp$BL = as.factor(cut(temp$d_N2,seq(0,1000000,block)))
  res<-t(unlist(lapply(lvl,function(x) err_D(temp[temp$BL %in% x,]))))
  ver<-c("H1","H2","H3","H123")
  vor<-c('exp','obs')
  vis<-paste(rep(vor, each = length(ver)), ver, sep = "")
  vars<-c(1:length(lvl))*block
  colnames(res)<-paste(rep(vars, each = length(vis)), vis, sep = "")
  # res<-as.data.frame(unlist(by(temp,temp$BL,function(x) err_D(x))))
  return(as.data.frame(res))
}

do_all<-function(file){
  ftree<-paste(file,'/spe_tree',sep='')
  tr<-read.tree(ftree)
  fD<-paste(file,'/data.txt',sep='')
  D<-read.table(fD,h=T,sep="\t")
  D$EV<-as.factor(paste(D$Donor,D$Recip, sep=""))
  D$FDR<-p.adjust(D$Pvalue,"fdr")
  fD3<-paste(file,'/data_D3.txt',sep='')
  D3<-read.table(fD3,h=T,sep="\t")
  D3$EV<-as.factor(paste(D3$Donor,D3$Recip, sep=""))
  res<-c(collesser(tr),err_D(D),err_D3(D3),block_err(D))#,tree_proba_err(tr))
  return(res)
}

adjus<-function(t) {
  return(p.adjust(t$Pvalue,"fdr"))
}

truc<-function(t) {
  res<-as.data.frame(rbind(by(t, t$Block, function(x) Ds_error_by_dist(x))))
  return(res)
}

Ds_error_by_dist<-function(temp){
  signi_D<-temp[temp$Adjust<=0.05,]
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
  alt<-nrow(signi_D[signi_D$EV == "N2P1"|signi_D$EV == "N2P2",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}

out<-function(temp){
  ### observed
  signi_D<-temp[temp$Adjust<=0.05,]
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
  alt<-nrow(signi_D[signi_D$EV == "OP1"|signi_D$EV == "OP2",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}

truc2<-function(t) {
  res<-as.data.frame(rbind(by(t, t$Block, function(x) out(x))))
  return(res)
}

out2<-function(temp){
  ### observed
  signi_D<-temp[temp$Adjust<=0.05,]
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
  alt<-nrow(signi_D[signi_D$EV == "P4P1"|signi_D$EV == "P4P2"|signi_D$EV == "P1P4"|signi_D$EV == "P2P4",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}

truc3<-function(t) {
  res<-as.data.frame(rbind(by(t, t$Block, function(x) out2(x))))
  return(res)
}


dist_relative<-function(x){
  N2<-x$dP1P4 - x$dP1P3
  res<-N2/max(N2)
  return(res)
}

tempdist_relative<-function(x){
  res<-x$dP1P4/ max(x$dP1P4)
  return(res)
}


# par(mfrow=c(1,3))

# error basic

dataDist<-read.table('new_05/data_D.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
dataDist$RD<-unlist(by(dataDist, dataDist$test, function(x) dist_relative(x)))
dataDist$Block = as.factor(cut(dataDist$RD,seq(0,1,0.1)))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))
ed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc(x)))
m_ed<-apply(ed,2,function(x) mean(na.omit(x)))

png(filename = "V2_distance_05.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(ed, main= "Amount of error by level of RDp3p4 (ext=0.5, block=0.1)", cex.axis=0.7, las = 1, varwidth = T,
  xlab="Levels of distance separating P3 and P4 (RDp3p4)", ylab="Proportion of mis-interpretation")
points(m_ed, col='red', pch = 18, cex = 2.3)
dev.off()

# error p4 and out
dataDist<-read.table('new_05/data_D.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
dataDist$RD<-unlist(by(dataDist, dataDist$test, function(x) dist_relative(x)))
dataDist$Block = as.factor(cut(dataDist$RD,seq(0,1,0.1)))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))
p4ed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc2(x)))
p4m_ed<-apply(p4ed,2,function(x) mean(na.omit(x)))
poed<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc3(x)))
pom_ed<-apply(poed,2,function(x) mean(na.omit(x)))


par(mfrow=c(1,3))
boxplot(ed, main= "Amount of error by level of RDp3p4 (ext=0.5, block=0.1)", las = 2 , varwidth = T)
points(m_ed, col='red', pch = 18, cex = 2.3)
boxplot(p4ed, main= "Amount of error by level of RDp3p4 (ext=0.5, block=0.1)", las = 2 , varwidth = T)
points(p4m_ed, col='blue', pch = 18, cex = 2.3)
boxplot(poed, main= "Amount of error by level of RDp3p4 (ext=0.5, block=0.1)", las = 2 , varwidth = T)
points(pom_ed, col='green', pch = 18, cex = 2.3)



# block=100000
# dataDist$Block = as.factor(cut(dataDist$dP1P4 - dataDist$dP1P3,seq(0,1000000,block)))
# edB<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc(x)))
# m_edB<-apply(ed,2,function(x) mean(na.omit(x)))
# boxplot(edB, main= "Amount of error by level of Dp1p2 (block=1e5, ext=0.9)", las = 2 , varwidth = T)
# points(m_edB, col='red', pch = 18, cex = 2)


ERR2<-function(temp){
  ### observed
  signi_D<-temp[temp$Adjust<=0.05,]
  signi_D<-temp
  true_pos<-nrow(signi_D[signi_D$EV == "P3P2"|signi_D$EV == "P2P3"|
    signi_D$EV == "P3P1"|signi_D$EV == "P1P3",])
    alt<-nrow(signi_D[signi_D$EV == "N2P1"|signi_D$EV == "N2P2",])
    # alt<-nrow(signi_D[signi_D$EV == "N2P1"|signi_D$EV == "N2P2"|signi_D$EV == "OP2"|signi_D$EV == "OP2"|signi_D$EV == "P4P2"|signi_D$EV == "P4P2",])
  obs_err<-alt/(true_pos+alt)
  return(obs_err)
}



dataDist<-read.table('data_D.txt',h=T)
time=read.table("TR_time",h=F)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))

p4ed<-by(dataDist, dataDist$test, function(x) ERR2(x))
err=as.data.frame(t(rbind(p4ed)))
plot(time[,2], err[,1])


# dist err D and D3

dataDist<-read.table('data_D.txt',h=T)
dataDist_D3<-read.table('data_D3.txt',h=T)


lvl<-levels(dataDist_D3$test)

D3_distancer<-function(x){
  temp<-dataDist[dataDist$test == as.character(unlist(x[1])),]
  p12<-temp[temp$P1 == as.character(unlist(x[2])) & temp$P2 == as.character(unlist(x[3])) | temp$P2 == as.character(unlist(x[2])) & temp$P3 == as.character(unlist(x[4])),][c(1:5),c("dP1P2","dP1P3")]
  true_p12<-as.numeric(names(sort(table(unlist(c(p12))),decreasing=TRUE))[1])
  p13<-temp[temp$P1 == as.character(unlist(x[2])) & temp$P3 == as.character(unlist(x[4])) | temp$P1 == as.character(unlist(x[2])) & temp$P4 == as.character(unlist(x[4])),][c(1:5),c("dP1P3","dP1P4")]
  true_p13<-as.numeric(names(sort(table(unlist(c(p13))),decreasing=TRUE))[1])
  res<-c("dP1P2" = true_p12, "dP1P3" = true_p13)
  return(res)
}


dataDist = dataDist[dataDist$test %in% c("test1001","test1002"),]
dataDist_D3 = dataDist_D3[dataDist_D3$test %in% c("test1001","test1002"),]
apply(dataDist_D3,1, function(x) D3_distancer(x))
