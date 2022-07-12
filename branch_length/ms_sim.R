#!/usr/bin/env Rscript
# Théo

args = commandArgs(trailingOnly=TRUE)

output="."
library(parallel)

library(reshape2)
library(ggplot2)

################################################################################



patterns_count<-function(sites){
  return(apply(sites, 2, function(x)
    if (all(x == c(1,0,0))){"aab"}
    else if (all(x == c(0,0,1))){"baa"}
    else if (all(x == c(0,1,0))){"aba"}
    else if (all(x == c(0,1,1))){"bba"}
    else{"other"}
  ))
}
#
#
# TT_count<-function(sim){
#   tt=read.tree(text=sim$trees[[1]][1])
#   dist<-cophenetic(tt)
#   tip=tt$tip.label
#   tip=tip[tip!="s1"]
#   dist=dist[tip,tip]
#   temp=apply(dist, 1 , sum)
#   out=names(temp[order(temp)][3])
#   if (out == "s2"){
#     topo="ABC"
#     pa=patterns_count(sim$seg_sites[[1]][c(2:4),])
#   }else if (out == "s3"){
#     topo="ACB"
#     site=sim$seg_sites[[1]][c(2:4),]
#     site=site[c(2,1,3),]
#     pa=patterns_count(site)
#   }else if (out == "s4"){
#     topo="BCA"
#     site=sim$seg_sites[[1]][c(2:4),]
#     site=site[c(3,1,2),]
#     pa=patterns_count(site)
#   }
#   T2=(1 / 1000) *  ((length(which(pa == 'aba')) + length(which(pa == 'baa')))/2)
#   T1=(1 / 1000) * (((length(which(pa == 'aba')) + length(which(pa == 'baa')))/2) + length(which(pa == 'bba')))
#   res=c("Topology"=topo, "T2"=T2, "T1"=T1)
#   return(res)
# }



TT_count<-function(sim){
  tt=read.tree(text=sim$trees[[1]][1])
  dist<-cophenetic(tt)
  temp=apply(dist, 1 , sum)
  out=names(temp[order(temp)][3])
  if (out == "s1"){
    topo="ABC"
    pa=patterns_count(sim$seg_sites[[1]])
  }else if (out == "s2"){
    topo="ACB"
    site=sim$seg_sites[[1]]
    site=site[c(2,1,3),]
    pa=patterns_count(site)
  }else if (out == "s3"){
    topo="BCA"
    site=sim$seg_sites[[1]]
    site=site[c(3,1,2),]
    pa=patterns_count(site)
  }
  T2=(1 / 1000) *  ((length(which(pa == 'aba')) + length(which(pa == 'baa')))/2)
  T1=(1 / 1000) * (((length(which(pa == 'aba')) + length(which(pa == 'baa')))/2) + length(which(pa == 'bba')))
  res=c("Topology"=topo, "T2"=T2, "T1"=T1)
  return(res)
}

std <- function(x) sd(x)/sqrt(length(x))

# # model <- coal_model(sample_size = c(1, 1, 1, 1), loci_number = 1, loci_length = 1, ploidy = 1) +
# # feat_mutation(1000, fixed_number = TRUE, locus_group = 'all') +
# # feat_pop_merge(1  , 4, 3) +
# # feat_pop_merge(1.2, 3, 2) +
# # feat_pop_merge(10  , 2, 1) +
# # # feat_migration(4000.0, pop_from = 3, pop_to = 2, symmetric = FALSE, time = 0, locus_group = 'all') +
# # # feat_migration(0, pop_from = 3, pop_to = 2, symmetric = FALSE, time = 0.000025, locus_group = 'all') +
# # sumstat_seg_sites() +
# # sumstat_trees()
#
# model <- coal_model(sample_size = c(1, 1, 1), loci_number = 1, loci_length = 1000, ploidy = 1) +
# # feat_mutation(1000, fixed_number = TRUE, locus_group = 'all') +
# feat_mutation(rate = 1000, model = 'IFS', fixed_number = FALSE, locus_group = 'all') +
# feat_pop_merge(0.5660125, 3, 2) +
# feat_pop_merge(1.44198, 2, 1) +
# feat_migration(0, symmetric = TRUE, time = 0, locus_group = 'all') +
# sumstat_seg_sites() +
# sumstat_trees()
# rep = simulate(model, nsim = 1000, cores = 8)
# TT<-as.data.frame(do.call('rbind',lapply(rep, function(x) TT_count(x))))
# TT$T1<-as.numeric(as.character(TT$T1))
# TT$T2<-as.numeric(as.character(TT$T2))
# res<-as.data.frame(do.call('rbind',by(TT, TT$Topology, function(x) {
#   res=c("mT1"=mean(x$T1), "mT2"=mean(x$T2), "seT1"=std(x$T1), "seT2"=std(x$T2))
#   return(res)
# })))
# par(mfrow=c(3,3))
# boxplot(TT$T1 ~ TT$Topology)
# boxplot(TT$T2 ~ TT$Topology)
# plot(res$mT1, ylim=c(min(res$mT2),max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n')
# axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
# points(res$mT2, pch = 19 , cex = 2 , col='blue')
# table(TT$Topology)
# legend("bottomright",
#   legend = c("T1", "T2"),
#   col = c("red", "blue"),
#   pch = c(19,19),
#   bty = "n",
#   pt.cex = 2,
#   cex = 1.2,
#   text.col = "black",
#   horiz = F ,
#   inset = c(0.1, 0.1))
#
# model <- coal_model(sample_size = c(1, 1, 1), loci_number = 1, loci_length = 1, ploidy = 1) +
# feat_mutation(rate = 1000, model = 'IFS', fixed_number = FALSE, locus_group = 'all') +
# # feat_mutation(1000, fixed_number = TRUE, locus_group = 'all') +
# feat_pop_merge(0.5660125, 3, 2) +
# feat_pop_merge(1.44198, 2, 1) +
# feat_migration(40000.0, pop_from = 3, pop_to = 1, symmetric = FALSE, time = 0, locus_group = 'all') +
# feat_migration(0, pop_from = 3, pop_to = 1, symmetric = FALSE, time = 0.0000025, locus_group = 'all') +
# feat_migration(40000.0, pop_from = 2, pop_to = 1, symmetric = FALSE, time = 0.01, locus_group = 'all') +
# feat_migration(0, pop_from = 2, pop_to = 1, symmetric = FALSE, time = 0.0100025, locus_group = 'all') +
# sumstat_seg_sites() +
# sumstat_trees()
# rep = simulate(model, nsim = 1000, cores = 8)
# TT<-as.data.frame(do.call('rbind',lapply(rep, function(x) TT_count(x))))
# TT$T1<-as.numeric(as.character(TT$T1))
# TT$T2<-as.numeric(as.character(TT$T2))
# res<-as.data.frame(do.call('rbind',by(TT, TT$Topology, function(x) {
#   res=c("mT1"=mean(x$T1), "mT2"=mean(x$T2), "seT1"=std(x$T1), "seT2"=std(x$T2))
#   return(res)
# })))
# boxplot(TT$T1 ~ TT$Topology)
# boxplot(TT$T2 ~ TT$Topology)
# plot(res$mT1, ylim=c(min(res$mT2),max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n')
# axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
# points(res$mT2, pch = 19 , cex = 2 , col='blue')
# table(TT$Topology)
# legend("bottomright",
#   legend = c("T1", "T2"),
#   col = c("red", "blue"),
#   pch = c(19,19),
#   bty = "n",
#   pt.cex = 2,
#   cex = 1.2,
#   text.col = "black",
#   horiz = F ,
#   inset = c(0.1, 0.1))
#
#
#


SEED=805
N_CORE=6
N_SIMULATION = 10
source(paste(output,"ms_command_sampled.R", sep ="/"))

rep = simulate(model, nsim = 1000, cores = 8)
TT<-as.data.frame(do.call('rbind',lapply(rep, function(x) TT_count(x))))
TT$T1<-as.numeric(as.character(TT$T1))
TT$T2<-as.numeric(as.character(TT$T2))
res<-as.data.frame(do.call('rbind',by(TT, TT$Topology, function(x) {
  res=c("mT1"=mean(x$T1), "mT2"=mean(x$T2), "seT1"=std(x$T1), "seT2"=std(x$T2))
  return(res)
})))
res
par(mfrow=c(2,3))
boxplot(TT$T1 ~ TT$Topology)
boxplot(TT$T2 ~ TT$Topology)
plot(res$mT1, ylim=c(min(res$mT2),max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n')
axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
points(res$mT2, pch = 19 , cex = 2 , col='blue')
legend("bottomright", legend = c("T1", "T2"), col = c("red", "blue"),
pch = c(19,19), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black",
horiz = F , inset = c(0.1, 0.1))
source(paste(output,"ms_command_ghost.R", sep ="/"))
rep = simulate(model, nsim = 1000, cores = 8)
TT<-as.data.frame(do.call('rbind',lapply(rep, function(x) TT_count(x))))
TT$T1<-as.numeric(as.character(TT$T1))
TT$T2<-as.numeric(as.character(TT$T2))
res<-as.data.frame(do.call('rbind',by(TT, TT$Topology, function(x) {
  res=c("mT1"=mean(x$T1), "mT2"=mean(x$T2), "seT1"=std(x$T1), "seT2"=std(x$T2))
  return(res)
})))
res
boxplot(TT$T1 ~ TT$Topology)
boxplot(TT$T2 ~ TT$Topology)
plot(res$mT1, ylim=c(min(res$mT2),max(res$mT1)), pch = 19 , cex = 2 , col='red', xaxt='n')
axis(side = 1, at = c(1:3), labels=c("ABC", "ACB", "BCA"), tick = F)
points(res$mT2, pch = 19 , cex = 2 , col='blue')
legend("bottomright", legend = c("T1", "T2"), col = c("red", "blue"),
pch = c(19,19), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black",
horiz = F , inset = c(0.1, 0.1))



#
#
#
# table(TT$Topology)
#
#
#
#
#
#
# par(mfrow=c(1,2))
# temp=TT[which(TT$Topology=="ABC"),]
#
# plot(sort(temp$T1) , ylim=c(min(temp[,c(2,3)]),max(temp[,c(2,3)])))
# points(sort(temp$T2))
# temp=TT[which(TT$Topology=="BCA"),]
#
# points(sort(temp$T1) , col='red')
# points(sort(temp$T2) , col='red')
#
#
#
#
# # par(mfrow=c(1,2))
#
#
#
#
# # lapply(rep, function(x) x$trees[[1]][1])
#
#
# aa=rep[[4]]
# aa$seg_sites[[1]]
# aa$trees[[1]]
# pa=aa$seg_sites[[1]]
#
# rep[[10]]$seg_sites[[1]]
# pa=patterns_count(rep[[10]]$seg_sites[[1]])
# length(which(pa == 'aba'))
# length(which(pa == 'baa'))
# length(which(pa == 'bba'))
#
# uniquer<-function(x){
#   if (ncol(x$seg_sites[[1]]) == 1) { # ajouter filtre sur site sum(seg) == 1
#     return(x$seg_sites[[1]][[1]])
#   }
# }
#
# sites<-do.call("cbind",mclapply(rep, function(x) uniquer(x),mc.cores = N_CORE))
#
# patterns<-apply(sites, 2, function(x)
#   if (sum(x) == 1){
#     if (x[1] == 1){"aab"}
#     else if (x[2] == 1){"aba"}
#     else if (x[3] == 1){"baa"}
#   }else if (sum(x) == 2){
#     if (x[2] == x[3]) {"bba"}
#     else {"other"}
#   }else {"other"}
# )
# table(patterns)
#
# T2 =
#
#
# # feat_migration(40000.0, pop_from = 3, pop_to = 2, symmetric = F, time = 0, locus_group = 'all') +
# # feat_migration(0, pop_from = 3, pop_to = 2, symmetric = F, time = 0.0000025, locus_group = 'all') +
# # feat_migration(1000, symmetric = T, time = 0, locus_group = 'all') +
# sumstat_trees() +
# source(paste(output,"ms_command_sampled.R", sep ="/"))
# rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
# single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))
# plot(single_trees[[30]])
# TTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
# TTs$T1<-as.numeric(as.character(TTs$T1))
# TTs$T2<-as.numeric(as.character(TTs$T2))
# par(mfrow=c(1,3))
# boxplot(TTs$T1 ~ TTs$Topology)
# boxplot(TTs$T2 ~ TTs$Topology)
# meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))
#
#
#
#
#
# which_topology<-function(tt){
#   dist<-cophenetic(tt)
#   if (tt$tip.label[1] == "s1"){
#     topo="ABC"
#     T1=dist["s3","s2"]/2
#     T2=dist["s3","s1"]/2
#   }else if (tt$tip.label[1] == "s3"){
#     topo="BCA"
#     T1=dist["s2","s1"]/2
#     T2=dist["s3","s1"]/2
#   }else if (tt$tip.label[1] == "s2"){
#     topo="ACB"
#     T1=dist["s1","s3"]/2
#     T2=dist["s2","s1"]/2
#   }
#   res=c("Topology"=topo, "T1"=T1, "T2"=T2)
#   return(res)
# }
# meaner<-function(df){
#   res=c("mT1"= mean(df[,2]), "mT2"= mean(df[,3]))
# }
#
# ################################################################################

SEED=805
N_CORE=8
N_SIMULATION = 1000
source(paste(output,"ms_command_ghost.R", sep ="/"))

rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))
cTTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
cTTs$T1<-as.numeric(as.character(cTTs$T1))
cTTs$T2<-as.numeric(as.character(cTTs$T2))
# boxplot(cTTs$T1 ~ cTTs$Topology, ylim = c(0, max(cTTs[,c(2,3)])))
# boxplot(cTTs$T2 ~ cTTs$Topology, ylim = c(0, max(cTTs[,c(2,3)])))
cmeanTTs<-as.data.frame(do.call('rbind',by(cTTs, cTTs$Topology, function(x) meaner(x))))
# plot(cmeanTTs$mT1, col="red", ylim=c(min(cmeanTTs),max(cmeanTTs)))
# points(cmeanTTs$mT2, col='blue')
write.table(cTTs,file="TT_ghost.txt",sep = "\t", row.names = F, append = F, quote=F)

cmeanTTs$Top<-row.names(cmeanTTs)
cmelTT<-melt(cmeanTTs)
cmelTT$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
  rgb(27, 158, 119, max=255))


par(mfrow=c(2,3))
# plot(melTT$value, col="black", pch = 23, bg=melTT$col, cex = 2,
#   xlim=c(0.8, 6.2),
#   frame.plot = FALSE)
#
# meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))
boxplot(cTTs$T1 ~ cTTs$Topology)
boxplot(cTTs$T2 ~ cTTs$Topology)

plot(cmeanTTs$mT1, xlab="", ylab="", col="red", xaxt='n', ylim=c(min(cmeanTTs$mT1),max(cmeanTTs$mT2)), cex = 3, pch = 19)
points(cmeanTTs$mT2, col='blue', cex = 3, pch = 19)
# table(TTs$Topology)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)
legend("bottomright",
  legend = c("T1", "T2"),
  col = c("red", "blue"),
  pch = c(19,19),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))




source(paste(output,"ms_command_sampled.R", sep ="/"))
rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))
TTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
TTs$T1<-as.numeric(as.character(TTs$T1))
TTs$T2<-as.numeric(as.character(TTs$T2))
# par(mfrow=c(1,3))
boxplot(TTs$T1 ~ TTs$Topology)
boxplot(TTs$T2 ~ TTs$Topology)
meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))

# png(filename = "Wi_I_T1_T2.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")

plot(meanTTs$mT1, xlab="", ylab="", col="red", xaxt='n', ylim=c(min(meanTTs$mT1),max(meanTTs$mT2)), cex = 3, pch = 19)
points(meanTTs$mT2, col='blue', cex = 3, pch = 19)
table(TTs$Topology)
# axis(side = 1, at = c(2,5), labels=c("T1", "T2"), tick = F, line = -2.5)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)
legend("bottomright",
  legend = c("T1", "T2"),
  col = c("red", "blue"),
  pch = c(19,19),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
# dev.off()











which_topology<-function(tt){
  dist<-cophenetic(tt)
  tip=tt$tip.label
  tip=tip[tip!="s1"]
  if (tip[1] == "s2"){
    topo="ABC"
    T1=dist["s3","s4"]/2

    T2=dist["s2","s3"]/2

  }else if (tip[1] == "s3"){
    topo="ACB"
    T1=dist["s2","s4"]/2

    T2=dist["s3","s2"]/2

  }else if (tip[1] == "s4"){
    topo="BCA"
    T1=dist["s2","s3"]/2

    T2=dist["s4","s3"]/2

  }
  res=c("Topology"=topo, "T1"=T1, "T2"=T2)
  return(res)
}
std <- function(x) sd(x)/sqrt(length(x))
meaner<-function(df){
  res=c("mT1"= mean(df[,2]), "mT2"= mean(df[,3]), "seT1" = std(df[,2]), "seT2" = std(df[,3]))
}

which_topology<-function(tt){
  dist<-cophenetic(tt)
  tip=tt$tip.label
  tip=tip[tip!="s1"]
  dist=dist[tip,tip]
  temp=apply(dist, 1 , sum)
  out=names(temp[order(temp)][3])
  if (out == "s2"){
    topo="ABC"
    T1=dist["s3","s4"]/2
    T2=dist["s2","s3"]/2
  }else if (out == "s3"){
    topo="ACB"
    T1=dist["s2","s4"]/2
    T2=dist["s3","s2"]/2
  }else if (out == "s4"){
    topo="BCA"
    T1=dist["s2","s3"]/2
    T2=dist["s4","s3"]/2
  }
  res=c("Topology"=topo, "T1"=T1, "T2"=T2)
  return(res)
}

SEED=805
N_CORE=8
N_SIMULATION = 1000
library('ape')
library('coala')
library('phyclust')
library('phangorn')
activate_ms(priority = 600)

model <- coal_model(sample_size = c(1, 1, 1, 1), loci_number = 1, loci_length = 1, ploidy = 1) +
feat_mutation(1, fixed_number = TRUE, locus_group = 'all') +
feat_pop_merge(1, 4, 3) +
feat_pop_merge(1.5, 3, 2) +
feat_pop_merge(5, 2, 1) +
feat_migration(40000.0, pop_from = 3, pop_to = 2, symmetric = F, time = 0, locus_group = 'all') +
feat_migration(0, pop_from = 3, pop_to = 2, symmetric = F, time = 0.0000025, locus_group = 'all') +
# feat_migration(1000, symmetric = T, time = 0, locus_group = 'all') +
sumstat_trees()
source(paste(output,"ms_command_sampled.R", sep ="/"))
rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))
plot(single_trees[[30]])
TTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
TTs$T1<-as.numeric(as.character(TTs$T1))
TTs$T2<-as.numeric(as.character(TTs$T2))
par(mfrow=c(1,3))
boxplot(TTs$T1 ~ TTs$Topology)
boxplot(TTs$T2 ~ TTs$Topology)
meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))

png(filename = "Wi_I_T1_T2.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")

plot(meanTTs$mT1, xlab="", ylab="", col="red", xaxt='n', ylim=c(min(meanTTs$mT1),max(meanTTs$mT2)), cex = 3, pch = 19)
points(meanTTs$mT2, col='blue', cex = 3, pch = 19)
table(TTs$Topology)
# axis(side = 1, at = c(2,5), labels=c("T1", "T2"), tick = F, line = -2.5)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)
legend("bottomright",
  legend = c("T1", "T2"),
  col = c("red", "blue"),
  pch = c(19,19),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
dev.off()


# write.table(TTs,file="TT_sampled.txt",sep = "\t", row.names = F, append = F, quote=F)
meanTTs$Top<-row.names(meanTTs)
melTT<-melt(meanTTs)
melTT$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
  rgb(27, 158, 119, max=255))

plot(melTT$value, col="black", pch = 23, bg=melTT$col, cex = 2,
  xlim=c(0.8, 6.2),
  frame.plot = FALSE)
axis(side = 1, at = c(2,5), labels=c("T1", "T2"), tick = F, line = -2.5)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)














library(reshape2)
library(ggplot2)

meaner<-function(df){
  res=c("mT1"= mean(df[,2]), "mT2"= mean(df[,3]))
}

ghost = read.table('TT_ghost.txt', h = T)
sampl = read.table('TT_sampled.txt', h = T)




png(filename = "Species_tree_T1_T2.png", width = 1200, height = 600, units = "px", pointsize = 12, bg = "white")

par(mfrow=c(1,2))

meanghostTTs<-as.data.frame(do.call('rbind',by(ghost, ghost$Topology, function(x) meaner(x))))
meanghostTTs$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
  rgb(27, 158, 119, max=255))
meanghostTTs$Top<-row.names(meanghostTTs)
meanghostTTs<-melt(meanghostTTs)
plot(meanghostTTs$value, col="black", pch = 23, bg=meanghostTTs$col, cex = 2,
  xlim=c(0.8, 6.2), ylim=c(0.8, max(meanghostTTs$value)),
  main = "With ghost species", xaxt='n')
abline(v=3.5)
axis(side = 1, at = c(2,5), labels=c("T1", "T2"), tick = F, line = -2.5)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)


meansamplTTs<-as.data.frame(do.call('rbind',by(sampl, sampl$Topology, function(x) meaner(x))))
meansamplTTs$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
  rgb(27, 158, 119, max=255))
meansamplTTs$Top<-row.names(meansamplTTs)
meansamplTTs<-melt(meansamplTTs)
plot(meansamplTTs$value, col="black", pch = 23, bg=meansamplTTs$col, cex = 2,
  xlim=c(0.8, 6.2), ylim=c(0.8, max(meanghostTTs$value)),
  main = "Only sampled species", xaxt='n')
abline(v=3.5)
axis(side = 1, at = c(2,5), labels=c("T1", "T2"), tick = F, line = -2.5)
axis(side = 1, at = c(1:6), labels=c("ABC", "ACB", "BCA", "ABC", "ACB", "BCA"), tick = F)
dev.off()



# # #
# # #
# # # meanTTs$Top<-row.names(meanTTs)
# # # melTT<-melt(meanTTs)
# # # melTT$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
# # #   rgb(27, 158, 119, max=255))
# # #
# # #
# # # plot(melTT$value, col="black", pch = 23, bg=melTT$col, cex = 2,
# # # # ylim=c(0, max(melTT$value)),
# # #   xlim=c(0.8, 6.2),
# # #   # bty="n",
# # #   frame.plot = FALSE,
# # #   # ylab = "", xlab = "", xaxt='n', yaxt='n'
# # #   )
# # # abline(h=0.875)
# # # abline(v=0.6)
# # #
# # # axis(side = 1, at = c(1,3.5, 6.5), labels=F)
# # # axis(side = 1, at = c(2.25,4.75), labels=c("T1", "T2"), tick = F)
# # # axis(side = 2)
# # #
# # # boxplot(x, frame.plot = FALSE,ylim=c(0,5))
# # # axis(side=1, pos=0, lwd.ticks=0)
# # # abline(h=0)
# # #
# # #
# # # ggplot(data = melTT, aes(x = variable, y = value, group = Top, fill = col)) +
# # #   geom_point()
# # #
# # #
# # #
# # #
# # #
# #
# #
# # # while(names(dev.cur()) !='null device') Sys.sleep(1)
# #
# # #
# # # tips<-tree$tip.label
# # # co <- rep("darkgrey", Nedge(tree))
# # # i <- which.edge(tree, tips[c(grep(c("uns"), tips),grep(c("ali"), tips))])
# # # co[i] <- "black"
# # # i <- which.edge(tree, tips[c(grep(c("ali"), tips))])
# # # co[i] <- "red"
# # # par(xpd = TRUE)
# # # plot(tree, edge.col = co, show.tip.label=F, direction="downwards", edge.width = 2)
# # #
# #
# # #
# # # t2<-keep.tip(tree, tips[c(grep(c("ali"), tips))])
# # #
# # # par(mfrow=c(1,2))
# # # plot(t2, edge.col = "red", show.tip.label=F, direction="downwards",
# # #   edge.width = 4, no.margin = F,
# # #   main = "Whithout concidering\nextinct lineages")
# # # plot(tree, edge.col = co, show.tip.label=F, direction="downwards",
# # #   edge.width = 2, no.margin = F,
# # #   main = "Concidering\nextinct lineages")
# # #
# #
# # #
# # # plot(tree, show.tip.label=F, direction="downwards", edge.width = 2)
