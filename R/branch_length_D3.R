#!/usr/bin/env Rscript
# Théo


require('ape')
require('phangorn')q

data<-read.table('data_D3.txt', h=T)
tree<-read.tree('spe_tree')
spnd<-c(tree$tip.label, tree$node.label)


how_many<-function(x){
  desc<-spnd[Descendants(tree,mrca.phylo(tree, c(as.numeric(as.character(x[1])),
    as.numeric(as.character(x[3])))), type="all")]
  return(c(
    "ali"=length(grep("ali>", desc)),
    "uns"=length(grep("uns>", desc)),
    "ext"=length(grep("ext>", desc))))
}

data[,c('ali','uns','ext')]=t(apply(data,1,function(x) how_many(x)))
data$R=3/(data$ali + data$uns + data$ext)

write.table(data,file="data_D3_v2.txt",sep = "\t", row.names = F, append = F, quote=F)


rm -rf test* spe_tree

for i in `seq 1301 1 1301`; do
   for j in 3; do
    sed "s/aaaa/$i/g" run_slurm.sh > temp
    sed -i "s/bbbb/$j/g" temp
    # sbatch temp
    bash temp
    rm temp
   done
done


sss test*/spe_tree
# cat test*/spe_tree > trees2
# sss trees2
#
#
# for i in test*;do
#   cd $i
#   /beegfs/data/soft/R-3.5.2/bin/Rscript ../scr.Rscript
#   cd ..
# done

column -t test2*/data_D3_v2.txt | grep "P[12]     P3\|P3     P[12]" | wc -l
column -t test2*/data_D3_v2.txt | grep "O      P[12]" | wc -l


rm data_D3.txt
files=(test*)
first=${files[0]}
files=("${files[@]:1}")
sed "s/^/`echo $first`\t/" $first/data_D3.txt > data_D3.txt
# sed "s/^/`echo $first`\t/" $first/data_Dfoil.txt > data_Dfoil.txt
name=`echo $first`
sed -i -e "0,/$name/ s/$name/test\t/" data_D3.txt
# sed -i -e "0,/$name/ s/$name/test\t/" data_Dfoil.txt
# process the remaining elements
for i in "${files[@]}"; do
  sed "1d" $i/data_D3.txt | sed "s/^/$i\t/" >> data_D3.txt
  # sed "1d" $i/data_Dfoil.txt | sed "s/^/$i\t/" >> data_Dfoil.txt
done
head data_D3.txt
# head data_Dfoil.txt

/beegfs/data/soft/R-3.5.2/bin/R

d<-read.table('data_D3.txt',h=T)
d$EV<-as.factor(paste(d$Donor,d$Recip, sep=""))
d<-d[which(abs(d$Z_score) >= 3),]

D3_error_by_param<-function(data){
  D3_H0<-length(which(data$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  D3_H1<-length(which(data$EV %in% c('OP1','OP2')))
  D3_err<-D3_H1/(D3_H0+D3_H1)
  return(c('err'=D3_err, "R"=nrow(data)))
}

a=by(d,as.factor(d$R), function(x) D3_error_by_param(x))
d$R2<-(d$ext+d$uns+d$ali)/d$allsp
step=0.05
xaxix=seq(0,1,step)
aa=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$R2<=x),]))))

# png(file = "D3_error_v3.png", width = 800, height = 800, res = 100)
svg(file = "D3_error_v3.svg", width = 8, height = 8)
plot(x=xaxix, y=aa$err, ylim=c(0,1),
  xlab= "",xaxt = "n",
  ylab= "", pch = 19, cex= 2, cex.lab= 1.3)
title(ylab="Proportion of erroneous interpretation of the D3",
  cex.lab=1.5, mgp=c(2.6,1,0))
title(xlab="Relative size of the ingroup tree (R)\ncomapred to the size of the complete tree",
  cex.lab=1.5, mgp=c(4,1,0))
axis(1, at=seq(0, 1, by=0.2), labels=c("0","<0.2","<0.4","<0.6","<0.8","<1"))
dev.off()









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







x=xaxix[12]
temp=d[which(d$R2<=x),]
temp[which(temp$EV %in% c('P3P1','P3P2','P1P3','P2P3')),]
temp[which(temp$EV %in% c('OP1','OP2')),]




D3_error_by_param2<-function(data){
  D3_H0<-length(which(data$EV %in% c('P3P1','P3P2','P1P3','P2P3')))
  D3_H1<-length(which(data$EV %in% c('OP1','OP2')))
  D3_err<-D3_H1/(D3_H0+D3_H1)
  # return(c('err'=D3_err, "R"=nrow(data)))
  return(D3_err)
}

test_drive<-function(data){
  return(as.data.frame(t(rbind(by(data,data$test, function(x) D3_error_by_param2(x)))), row.names=""))
}

d$R2<-3/(d$ext+d$uns+d$ali)

step=0.05
xaxix=seq(0,1,step)
temp = as.data.frame(do.call('cbind',lapply(xaxix,function(x) test_drive(d[which(d$R2>=x&d$R2<(x+step)),]))), col.names=xaxix)
colnames(temp) = xaxix
boxplot(temp)



xaxix=seq(0,0.99,0.05)
# d[which(d$R2>=x & d$R2 < (x+0.1)),]
aa=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$dN1>=x & d$dN1<x+0.05),]))))
plot(x=xaxix, y=aa$err, ylim=c(0,1))

aa=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$dN1>=x),]))))
plot(x=xaxix, y=aa$err, ylim=c(0,1))

d$ghost=d$uns+d$ext+d$ali
aa=as.data.frame(do.call('rbind',by(d,d$ghost, function(x) D3_error_by_param(x))))
plot(aa)


sort(unique(d$ghost))
aa=as.data.frame(do.call('rbind',lapply(sort(unique(d$ghost)),function(x) D3_error_by_param(d[which(d$ghost>=x),]))))
plot(y= aa$err, x = sort(unique(d$ghost)))




aa=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$R>=x),]))))
bb=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$dN1>=x),]))))
cc=as.data.frame(do.call('rbind',lapply(xaxix,function(x) D3_error_by_param(d[which(d$R2>=x),]))))

plot(x=xaxix, y=aa$err, col = "red")
points(x=xaxix, y=bb$err, col = "blue")
points(x=xaxix, y=cc$err, col = "green")




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

y=seq(0,0.99,0.1)

aa=as.data.frame(do.call('rbind',lapply(y,function(x) err_D3(d[which(d$RR>=x),]))))
