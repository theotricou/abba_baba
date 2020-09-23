############################################################################


library("ape")
library("phangorn")

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
    else
        return(NULL)
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
  signi_D<-temp
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

dist_relative<-function(x){
  N2<-x$dP1P4 - x$dP1P3
  res<-N2/x$dP1P4
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
  stemp<-temp[temp$Adjust >= 0.05,]
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



# plot
dataDist<-read.table('data_D.txt',h=T)
dataDist$EV<-as.factor(paste(dataDist$Donor,dataDist$Recip, sep=""))
dataDist$Adjust<-unlist(by(dataDist, dataDist$test, function(x) adjus(x)))
dataDist$N2_H.H<-as.factor(dataDist$N2_H.H)

# Error by test

res<-as.data.frame(do.call('rbind',by(dataDist, dataDist$test, function(x) err_D(x))))

write.table(as.data.frame(res$Derr_obs), "err_obs_09", sep = "\t", row.names = F, append = F, quote=F)





# Error by number of hidden lineage
res<-by(dataDist, dataDist$N2_H.H, function(x) Ds_error_by_dist(x))
a=as.data.frame(t(rbind(res)))
a$n<-as.numeric(unlist(rownames(a)))
r2<-summary(lm(a$res~a$n))$r.squared

# Error by distance N2
dataDist$RD<-unlist(by(dataDist, dataDist$test, function(x) dist_relative(x)))
dataDist$Block = as.factor(cut(dataDist$RD,seq(0,1,0.1)))
red<-do.call(rbind,by(dataDist, dataDist$test, function(x) truc(x)))
m_red<-apply(red,2,function(x) mean(na.omit(x)))


png(filename = "Distance_00.png", width = 1200, height = 600, units = "px", pointsize = 12, bg = "white")
par(mfrow=c(1,2))
plot(res ,main=paste('Observed erro by hidden diversity (r2=', round(r2,5),')',sep=""),
  xlab="Number of hidden lineage branching on N2", ylab="Proportion of observed error", pch=19)
abline(lm(a$res~a$n),col='red')
boxplot(red, main= "Obverved error by the distance N2", cex.axis=0.7, las = 1, varwidth = T,
  xlab="Levels of distance N2", ylab="Obverved error")
points(m_red, col='red', pch = 18, cex = 2)
dev.off()


# colles plot

d0<-read.table("/beegfs/data/tricou/abba_baba/proba_00")
d5<-read.table("/beegfs/data/tricou/abba_baba/proba_05")
d9<-read.table("/beegfs/data/tricou/abba_baba/proba_09")

maxcol_c<-max(d0$V2,d5$V2,d9$V2)
maxcol_e<-max(d0$V3,d5$V3,d9$V3)
mincol_e<-min(d0$V3,d5$V3,d9$V3)
maxer_e<-max(d0$V1,d5$V1,d9$V1)
miner_e<-min(d0$V1,d5$V1,d9$V1)

png(filename = "Colless_00.png", width = 1000, height = 600, units = "px", pointsize = 12, bg = "white")
par(mfrow=c(1,1))
plot(d0$V3,d0$V1, xlim=c(mincol_e,maxcol_e), ylim=c(miner_e,maxer_e), col="blue", pch = 19, cex=0.5,
  main="Probability of error by extant tree imbalance (colless index)",
  xlab="Colless index", ylab="Probability of error")
legend("topleft",legend=c("ext=0.0", "ext=0.5","ext=0.9"),
       col=c("blue", "green", "red"), lty=1, cex=2)
points(d5$V3,d5$V1, xlim=c(mincol_e,maxcol_e), ylim=c(miner_e,maxer_e), col="green", pch = 19, cex=0.5)
points(d9$V3,d9$V1, xlim=c(mincol_e,maxcol_e), ylim=c(miner_e,maxer_e), col="red", pch = 19, cex=0.5)
abline(lm(d0$V1~d0$V3),col='blue')
abline(lm(d5$V1~d5$V3),col='green')
abline(lm(d9$V1~d9$V3),col='red')
dev.off()



######### test time vs error


par(mfrow=c(1,1))
err0<-read.table('multi_tree_0/err_obs_00',h=T)
tim0<-read.table('time_00',h=F)
d0=na.omit(cbind(tim0,err0))


err5<-read.table('multi_tree_5/err_obs_05',h=T)
tim5<-read.table('time_05',h=F)
d5=na.omit(cbind(tim5,err5))

err9<-read.table('multi_tree_9/err_obs_09',h=T)
tim9<-read.table('time_09',h=F)
d9=na.omit(cbind(tim9,err9))


plot(d0, col='green', ylim =c(min(c(d0[,2],d5[,2],d9[,2])),max(c(d0[,2],d5[,2],d9[,2]))))
points(d5, col='blue')
points(d9, col='red')
