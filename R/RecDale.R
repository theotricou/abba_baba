#!/usr/bin/env Rscript
# Théo
# RecDale creation des quartet



quarteter<-function(x){
  if (x[16] %in% c("OP1","N2P1","P1P4","P4P1","P3P2","P2P3")){ # for D>0 ==> ABBA
    return(paste("(",x[1],",",x[4],",(",x[2],",",x[3],"));", sep=""))
  } else if (x[16] %in% c("OP2","N2P2","P2P4","P4P2","P3P1","P1P3")){ # for D<0 ==> BABA)
    return(paste("(",x[2],",",x[4],",(",x[1],",",x[3],"));", sep=""))
  } else{
    return(paste("(",x[4],",",x[3],",(",x[1],",",x[2],"));", sep=""))
  }
}


QMCqarteter<-function(x){
  if (x[16] %in% c("OP1","N2P1","P1P4","P4P1","P3P2","P2P3")){ # for D>0 ==> ABBA
    return(paste(x[4],",",x[1],"|",x[3],",",x[2], sep=""))
  } else if (x[16] %in% c("OP2","N2P2","P2P4","P4P2","P3P1","P1P3")){ # for D<0 ==> BABA)
    return(paste(x[4],",",x[2],"|",x[1],",",x[3], sep=""))
  } else{
    return(paste(x[4],",",x[3],"|",x[1],",",x[2], sep=""))
  }
}


paste(x[4],",",x[1],"|",x[3],",",x[2], sep="")
paste(x[4],",",x[2],"|",x[1],",",x[3], sep="")
paste(x[4],",",x[3],"|",x[1],",",x[2], sep="")


dataDist<-read.table('data.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))




trees=as.data.frame(apply(dataDist,1,function(x) quarteter(x)))
trees=as.data.frame(apply(dataDist,1,function(x) QMCqarteter(x)))

write.table(trees, "quartet_trees", sep = "\t", row.names = F,col.names = F, append = F, quote=F)


dataDist<-read.table('data.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
