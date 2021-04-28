#!/usr/bin/env Rscript
# Théo

args = commandArgs(trailingOnly=TRUE)

output="."
source(paste(output,"ms_command_ghost.R", sep ="/"))
library(parallel)

################################################################################

which_topology<-function(tt){
  dist<-cophenetic(tt)
  if (tt$tip.label[1] == "s1"){
    topo="ABC"
    T1=dist["s3","s2"]/2
    T2=dist["s3","s1"]/2
  }else if (tt$tip.label[1] == "s3"){
    topo="BCA"
    T1=dist["s2","s1"]/2
    T2=dist["s3","s1"]/2
  }else if (tt$tip.label[1] == "s2"){
    topo="ACB"
    T1=dist["s1","s3"]/2
    T2=dist["s2","s1"]/2
  }
  res=c("Topology"=topo, "T1"=T1, "T2"=T2)
  return(res)
}
meaner<-function(df){
  res=c("mT1"= mean(df[,2]), "mT2"= mean(df[,3]))
}

################################################################################

SEED=805
N_CORE=8
N_SIMULATION = 3000

# X11()


rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)

single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))


TTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
TTs$T1<-as.numeric(as.character(TTs$T1))
TTs$T2<-as.numeric(as.character(TTs$T2))
par(mfrow=c(2,3))
boxplot(TTs$T1 ~ TTs$Topology)
boxplot(TTs$T2 ~ TTs$Topology)




meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))
plot(meanTTs$mT1, col="red", ylim=c(min(meanTTs),max(meanTTs)))
points(meanTTs$mT2, col='blue')
write.table(TTs,file="TT_ghost.txt",sep = "\t", row.names = F, append = F, quote=F)


source(paste(output,"ms_command_sampled.R", sep ="/"))

rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)

single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))

TTs<-as.data.frame(do.call('rbind',lapply(single_trees, function(x) which_topology(x))))
TTs$T1<-as.numeric(as.character(TTs$T1))
TTs$T2<-as.numeric(as.character(TTs$T2))
boxplot(TTs$T1 ~ TTs$Topology)
boxplot(TTs$T2 ~ TTs$Topology)

meanTTs<-as.data.frame(do.call('rbind',by(TTs, TTs$Topology, function(x) meaner(x))))
plot(meanTTs$mT1, col="red", ylim=c(min(meanTTs),max(meanTTs)))
points(meanTTs$mT2, col='blue')


write.table(TTs,file="TT_sampled.txt",sep = "\t", row.names = F, append = F, quote=F)




library(reshape2)
library(ggplot2)


meanTTs$Top<-row.names(meanTTs)
melTT<-melt(meanTTs)
melTT$col=c(rgb(217, 94, 2, max=255), rgb(117, 112, 179, max=255),
  rgb(27, 158, 119, max=255))


plot(melTT$value, col="black", pch = 23, bg=melTT$col, cex = 2,
# ylim=c(0, max(melTT$value)),
  xlim=c(0.8, 6.2),
  # bty="n",
  frame.plot = FALSE,
  ylab = "", xlab = "", xaxt='n', yaxt='n'
  )
abline(h=0.875)
abline(v=0.6)

axis(side = 1, at = c(1,3.5, 6.5), labels=F)
axis(side = 1, at = c(2.25,4.75), labels=c("T1", "T2"), tick = F)
axis(side = 2)

boxplot(x, frame.plot = FALSE,ylim=c(0,5))
axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)


ggplot(data = melTT, aes(x = variable, y = value, group = Top, fill = col)) +
  geom_point()







# while(names(dev.cur()) !='null device') Sys.sleep(1)

#
# tips<-tree$tip.label
# co <- rep("darkgrey", Nedge(tree))
# i <- which.edge(tree, tips[c(grep(c("uns"), tips),grep(c("ali"), tips))])
# co[i] <- "black"
# i <- which.edge(tree, tips[c(grep(c("ali"), tips))])
# co[i] <- "red"
# par(xpd = TRUE)
# plot(tree, edge.col = co, show.tip.label=F, direction="downwards", edge.width = 2)
#

#
# t2<-keep.tip(tree, tips[c(grep(c("ali"), tips))])
#
# par(mfrow=c(1,2))
# plot(t2, edge.col = "red", show.tip.label=F, direction="downwards",
#   edge.width = 4, no.margin = F,
#   main = "Whithout concidering\nextinct lineages")
# plot(tree, edge.col = co, show.tip.label=F, direction="downwards",
#   edge.width = 2, no.margin = F,
#   main = "Concidering\nextinct lineages")
#

#
# plot(tree, show.tip.label=F, direction="downwards", edge.width = 2)
