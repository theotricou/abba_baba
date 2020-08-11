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

  # res<-c("Derr_exp"=D_err_exp_H1,"Derr_obs"=sD_err_exp_H1)

  res<-c("Derr_exp"=D_err_exp_H1,"DH2err_exp"=D_err_exp_H2,"DH3err_exp"=D_err_exp_H3,"DH123err_exp"=D_err_exp_H123,
    "Derr_obs"=sD_err_exp_H1,"DH2err_obs"=sD_err_exp_H2,"DH3err_obs"=sD_err_exp_H3,"DH123err_obs"=sD_err_exp_H123
  )

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
  # created by Damien
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  }
}
getquatuors<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}
#TEST IF TWO INTERVALS OVERLAP
overlap<-function(ab,cd) {
    #true if the two intervals overlap
    return(!(ab[2]<cd[1] | cd[2]<ab[1]))
}
getstartstopofedges<-function(mat,xx) {
    if (!is.null(mat)) {
        return(t(apply(mat,1,function(x,xx) c(x,xx[x]),xx=xx)))
    }
    else return(NULL)
}

plottreeandinfo<-function(tr, quat, returnxy=TRUE, plot=TRUE) {
    plot(tr, plot=plot)
    xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
    # p1p2<-getMRCA(tr,c(quat[1],quat[2]))
    # p1p3<-getMRCA(tr,c(quat[1],quat[3]))
    # p1p4<-getMRCA(tr,c(quat[1],quat[4]))
    # points(xy$xx[c(p1p2, p1p3,p1p4)], xy$yy[c(p1p2, p1p3,p1p4)], col=c("red","blue","green"), pch=19)
    # abline(v=xy$xx[c(p1p2, p1p3,p1p4)], col=c("red","blue","green"))
    #
    # p1nb<-which(tr$tip.label==quat[1])
    # p2nb<-which(tr$tip.label==quat[2])
    # p3nb<-which(tr$tip.label==quat[3])
    # p4nb<-which(tr$tip.label==quat[4])
    # segments(xy$xx[p1nb], xy$yy[p1nb], xy$xx[p1p2], xy$yy[p1p2], col="green")
    # segments(xy$xx[p2nb], xy$yy[p2nb], xy$xx[p1p2], xy$yy[p1p2], col="green")
    # segments(xy$xx[p1p2], xy$yy[p1p2], xy$xx[p1p3], xy$yy[p1p3], col="green")
    # segments(xy$xx[p3nb], xy$yy[p3nb], xy$xx[p1p3], xy$yy[p1p3], col="green")
    # segments(xy$xx[p1p3], xy$yy[p1p3], xy$xx[p1p4], xy$yy[p1p4], col="green")
    # segments(xy$xx[p4nb], xy$yy[p4nb], xy$xx[p1p4], xy$yy[p1p4], col="green")
    if (returnxy) return(xy)
}

probas.donor.receiver<-function(quat, tr, xy) { #xy are x and y coordinates in the plot
#    plottreeandinfo(tr, quat)
    p1<-quat[1]
    p2<-quat[2]
    p3<-quat[3]
    p4<-quat[4]
    #################
    ### P1 AND P2 ###
    #################
    p1p2<-mrca.phylo(tr, c(p1,p2))
    p1nb<-which(tr$tip.label==p1)
    p2nb<-which(tr$tip.label==p2)
    desc<-sapply(Children(tr, p1p2), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p2nb, desc[[1]])) desc<-rev(desc) ##so that p1 remains in first position.
    bl1.donor<-sum(tr$edge.length[match(desc[[1]],tr$edge[,2])])
    bl2.donor<-sum(tr$edge.length[match(desc[[2]],tr$edge[,2])])
    ##WE STORE ALL POSSIBLE BRANCHES
    P1DONOR<-tr$edge[match(desc[[1]],tr$edge[,2]),, drop=FALSE]
    P2DONOR<-tr$edge[match(desc[[2]],tr$edge[,2]),,drop=FALSE]

    ## NOW WE LOOK AT THE SUM OF BL AS RECEIVERS. fOR AN INTROGRESSION TO BE DETECTED,
    ## THE RECEIVING BRANCH MUST BE AN ANCESTOR OF P1 (RESP. P2) BUT NO OLDER THAN THE MRCA OF P1P2
    p1.up<-c(p1nb, Ancestors(tr, p1nb))
    p1.up.small<-p1.up[1:(which(p1.up==p1p2)-1)]
    p2.up<-c(p2nb, Ancestors(tr, which(tr$tip.label==p2)))
    p2.up.small<-p2.up[1:(which(p2.up==p1p2)-1)]
    bl1.receiver = sum(tr$edge.length[match(p1.up.small,tr$edge[,2])])
    bl2.receiver = sum(tr$edge.length[match(p2.up.small,tr$edge[,2])])
    ##WE STORE ALL POSSIBLE BRANCHES
    P1RECEIVER<-tr$edge[match(p1.up.small,tr$edge[,2]),, drop=FALSE]
    P2RECEIVER<-tr$edge[match(p2.up.small,tr$edge[,2]),, drop=FALSE]

    #################
    ###     P3    ###
    #################
    ## FOR P3 and P4 IT IS A BIT MORE TRICKY BECAUSE WE MUST EXCLUDE BRANCHES THAT END BEFORE MRCA(P1,P2)
    ## AND INTEGRATE PART OF THOSE THAT ARE CROSSING THIS LINE
    ## NOTE: P4 WILL NEVER BE THE RECEIVER BECAUSE IT IS GHOST.
    p1p3<-mrca.phylo(tr, c(p1,p3)) #same than p2p3
    p3nb<-which(tr$tip.label==p3)
    desc1.3<-sapply(Children(tr, p1p3), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p1nb, desc1.3[[1]])) desc1.3<-rev(desc1.3)
    desc1.3<-desc1.3[[1]]
    subedge<-tr$edge[match(desc1.3,tr$edge[,2]),, drop=FALSE]
    ##
    xp1p2<-xy$xx[p1p2] #time of p1p2 speciation. p3 donor or receiver must be more recent.
    xp1<-xy$xx[p1nb]
    ##
    subedgetime<-cbind(subedge, t(apply(subedge,1,function(x,xx) xx[x], xx=xy$xx)))
    subedgetime<-subedgetime[(subedgetime[,3]>xp1p2)|(subedgetime[,4]>xp1p2),, drop=FALSE] ##select those that have one side at least in the range encounterd by p1p2
    subedgetime[subedgetime[,3]<xp1p2,3]<-xp1p2 #change start time of the edges traversing the line to only count possible bl.
    bl3.donor<-sum(subedgetime[,4]-subedgetime[,3])
    ##WE STORE ALL POSSIBLE BRANCHES
    P3DONOR<-subedge

    p3.up<-c(p3nb, Ancestors(tr, p3nb))
    p3.up.subedgetime<-match(p3.up, subedgetime[,2])
    subedgetime.up<-subedgetime[p3.up.subedgetime[!is.na(p3.up.subedgetime)],, drop=FALSE]
    bl3.receiver<-sum(subedgetime.up[,4]-subedgetime.up[,3])
    ##WE STORE ALL POSSIBLE BRANCHES
    P3RECEIVER<-subedgetime.up[,1:2, drop=FALSE]
    #################
    ###     P4    ###
    #################
    ## P4 represents a possible ghost lineage that branches between the ancestor of (P1P2P3) and the ancestor of (P1P2P3P4)
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

    ###############################################################################
    ### STORE DISTANCE BETWEEN OUTGROUP AND P1P2P3 (in proportion of full size) ###
    ###############################################################################
    relativedist2outgroup<-(xy$xx[p1p3]-xy$xx[p1p4])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top3<-(xy$xx[p1nb]-xy$xx[p1p3])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top2<-(xy$xx[p1nb]-xy$xx[p1p2])/(xy$xx[p1nb]-xy$xx[p1p4])
    ###############################
    ###  COMPUTE PROBABILITIES  ###
    ###############################
    treesize<-sum(tr$edge.length)
    pP1d<-bl1.donor/treesize
    pP2d<-bl2.donor/treesize
    pP3d<-bl3.donor/treesize
    pP4d<-bl4.donor/treesize
    pP1r<-bl1.receiver/treesize
    pP2r<-bl2.receiver/treesize
    pP3r<-bl3.receiver/treesize

    ############################################
    ### COMPUTE EDGE-EDGE POSSIBLE TRANSFERS ###
    ############################################
    DONOR.RECEIVER<-list(P1DONOR=P1DONOR, P2DONOR=P2DONOR,P3DONOR=P3DONOR,P4DONOR=P4DONOR,P1RECEIVER=P1RECEIVER,P2RECEIVER=P2RECEIVER,P3RECEIVER=P3RECEIVER)
    DONOR.RECEIVER<-lapply(DONOR.RECEIVER, getstartstopofedges, xx=xy$xx)
    ##THEN FOR EACH POSSIBLE TRABSFERSTORY (P1P3, P3P1, etc.) we only keep those possible in terms of temporal co-occurence of edges.
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
    # P is the list of probas computed from probas.donor.receiver() function.
    #The proba of detecting a false positiveis the proba that P1 is donor AND P3 is receiver, or P2 is donor and P3 is receiver, or...
    proba.true.positive<-(P$pP1d*P$pP3r)+(P$pP3d*P$pP1r)+(P$pP2d*P$pP3r)+(P$pP3d*P$pP2r)
    proba.false.positive<-(P$pP4d*P$pP2r)+(P$pP4d*P$pP1r)
    return(list(tp=proba.true.positive, fp=proba.false.positive))
}

tree_proba_err<-function(tr){
  extantsampled<-tr$tip.label[grep("l>",tr$tip.label)]
  extant.tr<-keep.tip(tr, extantsampled)
  quatuors <- getquatuors(extant.tr)

  xy<-plottreeandinfo(tr, quatuors[1,], returnxy=TRUE, plot=T)

  #get all probas for all quaturos
  PROBAS<-apply(quatuors, 1, probas.donor.receiver, tr=tr, xy=xy)

  ###POST-PROCESSING
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

dir<-c('one_tree/','one_tree5/','one_tree9/')
dir<-c('sim_0/','sim_5/','sim_9/')
# dir<-list.files(path = ".",pattern = 'one_tree')
# dir<-c('sim_0/')
#
#
# RES0<-as.data.frame(do.call(rbind,lapply(list.files(path = dir[1],
#    pattern = 'test'),function(x) unlist(do_all(paste(dir[1],x,sep="")))))
# )
#
# RES5<-as.data.frame(do.call(rbind,lapply(list.files(path = dir[2],
#    pattern = 'test'),function(x) unlist(do_all(paste(dir[2],x,sep="")))))
# )
# RES9<-as.data.frame(do.call(rbind,lapply(list.files(path = dir[3],
#    pattern = 'test'),function(x) unlist(do_all(paste(dir[3],x,sep="")))))
# )

#
# RES=RES0
# RES=RES5
# RES=RES9


RES<-as.data.frame(do.call(rbind,mclapply(list.files(path=dir[1],pattern = 'test'),
function(x) unlist(do_all(paste(dir[1],x,sep=""))), mc.cores=20)))


pdf(file = "D_D3_ERROR.pdf", width = 16, height = 8)

par(mfrow=c(1,2))
plot(RES$Derr_exp+1, RES$D3err_exp, main='Expected D error vs D3 error',
xlab='D error', ylab="D3 erro")
plot(RES$Derr_obs, RES$D3err_obs, main='Observed D error vs D3 error',
xlab='D error', ylab="D3 erro")
dev.off()



write.table(RES,file="RES_sim_00.txt", sep = "\t", row.names = F, append = F, quote=F)

RES<-as.data.frame(do.call(rbind,mclapply(list.files(path=dir[2],pattern = 'test'),
function(x) unlist(do_all(paste(dir[2],x,sep=""))), mc.cores=20)))

write.table(RES,file="RES_sim_05.txt", sep = "\t", row.names = F, append = F, quote=F)

RES<-as.data.frame(do.call(rbind,mclapply(list.files(path=dir[3],pattern = 'test'),
function(x) unlist(do_all(paste(dir[3],x,sep=""))), mc.cores=20)))

write.table(RES,file="RES_sim_09.txt", sep = "\t", row.names = F, append = F, quote=F)




cbind.fill <- function(...){
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x)
        rbind(x, matrix(, n-nrow(x), ncol(x)))))
}




### D err
pdf(file = "ONE_D_ERR.pdf", width = 10, height = 8)
d0<-read.table('RES_one_00.txt', h=T)
d5<-read.table('RES_one_05.txt', h=T)
d9<-read.table('RES_one_09.txt', h=T)
data<-as.data.frame(cbind.fill(d0$Derr_exp,d0$Derr_obs,d5$Derr_exp,d5$Derr_obs,d9$Derr_exp,d9$Derr_obs))
boxplot(data,xaxt="n",
  ylab="D error rate", xlab="Extinction rate",
  main="D statistic error rate, expected and observed, by extinction rate"
)
points(apply(data, 2, function(x) mean(na.omit(x))),col='red',pch=18, cex = 2)
points(c(d0$proba[1],d0$proba[1],d5$proba[1],d5$proba[1],d9$proba[1],d9$proba[1]),col='blue',pch=18, cex = 2)
vis<-c('exp\n','obs\n')
vars<-c('ext=0','ext=.5','ext=.9')
labs<-paste(vis, rep(vars, each = length(vis)), sep = "")
axis(side = 1, pos=-.05, lty=0,at = seq(1:6),labels = labs, cex.axis = 1)
dev.off()




### D3 err
pdf(file = "ONE_D3_ERR.pdf", width = 10, height = 8)
d0<-read.table('RES_one_00.txt', h=T)
d5<-read.table('RES_one_05.txt', h=T)
d9<-read.table('RES_one_09.txt', h=T)
data<-as.data.frame(cbind.fill(d0$D3err_exp,d0$D3err_obs,d5$D3err_exp,d5$D3err_obs,d9$D3err_exp,d9$D3err_obs))
boxplot(data,xaxt="n",
  ylab="D3 error rate", xlab="Extinction rate",
  main="D3 statistic error rate, expected and observed, by extinction rate"
)
points(apply(data, 2, function(x) mean(na.omit(x))),col='red',pch=18, cex = 2)
# points(c(d0$proba[1],d0$proba[1],d5$proba[1],d5$proba[1],d9$proba[1],d9$proba[1]),col='blue',pch=18, cex = 2)
vis<-c('exp\n','obs\n')
vars<-c('ext=0','ext=.5','ext=.9')
labs<-paste(vis, rep(vars, each = length(vis)), sep = "")
axis(side = 1, pos=-.05, lty=0,at = seq(1:6),labels = labs, cex.axis = 1)
dev.off()







# Error by Quat, for single tree simulation

do_quat<-function(file){
  fD<-paste(file,'/data.txt',sep='')
  D<-read.table(fD,h=T,sep="\t")
  D$EV<-paste(D$Donor,D$Recip, sep="")

  # fD3<-paste(file,'/data_D3.txt',sep='')
  # D3<-read.table(fD3,h=T,sep="\t")
  # D3$EV<-as.factor(paste(D3$Donor,D3$Recip, sep=""))

  return(D$EV)
}
Err_quat<-function(temp){
  D_H0<-length(which(temp %in% c('P3P1','P3P2','P1P3','P2P3')))
  D_H1<-length(which(temp %in% c('N2P1','N2P2')))
  D_H2<-length(which(temp %in% c('OP1','OP2')))
  D_H3<-length(which(temp %in% c('P4P1','P4P2')))

  D_err_exp_H1<-D_H1/(D_H0+D_H1)
  D_err_exp_H2<-D_H2/(D_H0+D_H2)
  D_err_exp_H3<-D_H3/(D_H0+D_H3)
  D_err_exp_H123<-(D_H1+D_H2+D_H3)/(D_H0+D_H1+D_H2+D_H3)

  res<-c("DH1err_exp"=D_err_exp_H1,"DH2err_exp"=D_err_exp_H2,
    "DH3err_exp"=D_err_exp_H3,"DH123err_exp"=D_err_exp_H123)
  return(res)
}



quat_proba_err<-function(tr){
  extantsampled<-tr$tip.label[grep("l>",tr$tip.label)]
  extant.tr<-keep.tip(tr, extantsampled)
  quatuors <- getquatuors(extant.tr)

  xy<-plottreeandinfo(tr, quatuors[1,], returnxy=TRUE, plot=T)

  #get all probas for all quaturos
  PROBAS<-apply(quatuors, 1, probas.donor.receiver, tr=tr, xy=xy)

  ###POST-PROCESSING
  PPP<-lapply(PROBAS, ProbaFalsePositive)
  RES<-do.call(rbind, lapply(PPP, unlist))
  result<-RES[,2]/(RES[,1]+RES[,2])
  # data<-c("proba"=mean(na.omit(result)))
  return(result)
}


dd<-dir[3]

RESt<-as.data.frame(do.call(rbind,mclapply(list.files(path=dd,pattern = 'test'),
  function(x) do_quat(paste(dd,x,sep="")), mc.cores=10))
)
ftree<-paste(dd,list.files(path=dd,pattern = 'test')[1],'/spe_tree',sep='')
tr<-read.tree(ftree)
truc<-cbind(quat_proba_err(tr),as.data.frame(do.call(rbind,lapply(RESt,function(x) Err_quat(x)))))


pdf(file = "QUAT_ERR_09.pdf", width = 8, height = 8)

r2=summary(lm(truc[,2]~0+truc[,1]))$r.squared
plot(truc[,c(1,2)],ylim=c(0,1),xlim=c(0,1),
main=paste('Quatuors D stat error rate (extinction rate = 0) (r2=',r2,')',sep=""),
xlab="Probability of error", ylab="Simulated error rate")
lines(x = c(0,1),y = c(0,1),col='blue')
abline(lm(truc[,2]~0+truc[,1]),col='red')

dev.off()




#
#
# pdf(file = "COLLESS_PROBA.pdf", width = 16, height = 8)
#
# par(mfrow=c(1,2))
# r2=summary(lm(d$proba~0+d$c_col))$r.squared
# plot(d$c_col, d$proba, main=paste('Complete tree colless vs D stat error probability (r2=',r2,')',sep=""),
# xlab='Colless score', ylab="Probability of error")
# abline(lm(d$proba~0+d$c_col),col='red')
#
# r2=summary(lm(d$proba~0+d$e_col))$r.squared
# plot(d$e_col, d$proba, main=paste('Extant tree colless vs D stat error probability (r2=',r2,')',sep=""),
# xlab='Colless score', ylab="Probability of error")
# abline(lm(d$proba~0+d$e_col),col='red')
#
# dev.off()
#
#


# library("plyr")
#
# # D error + extinction
# pdf(file = "ONE_D_ERR.pdf", width = 10, height = 8)
#
# data=t(rbind.fill(as.data.frame(t(RES0[,c(3,4)])), as.data.frame(t(RES5[,c(3,4)])),
#   as.data.frame(t(RES9[,c(3,4)]))))
# boxplot(data,xaxt="n",
#   ylab="D error rate", xlab="Extinction rate",
#   main="D statistic error rate, expected and observed, by extinction rate")
#
# points(apply(data, 2, function(x) mean(na.omit(x))),col='red',pch=18, cex = 2)
# vis<-c('exp\n','obs\n')
# vars<-c('ext_00','ext_05','ext_09')
# labs<-paste(vis, rep(vars, each = length(vis)), sep = "")
# axis(side = 1, pos=-.05, lty=0,at = seq(1:6),labels = labs, cex.axis = 1)
#
# dev.off()
#
#
#
#
# # D3 error + extinction
# pdf(file = "ONE_D3_ERR.pdf", width = 10, height = 8)
#
# data=t(rbind.fill(as.data.frame(t(RES0[,c(5,6)])), as.data.frame(t(RES5[,c(5,6)])),
#   as.data.frame(t(RES9[,c(5,6)]))))
# boxplot(data,xaxt="n",
#   ylab="D3 error rate", xlab="Extinction rate",
#   main="D3 error rate, expected and observed, by extinction rate")
#
# points(apply(data, 2, function(x) mean(na.omit(x))),col='red',pch=18, cex = 2)
# vis<-c('exp\n','obs\n')
# vars<-c('ext_00','ext_05','ext_09')
# labs<-paste(vis, rep(vars, each = length(vis)), sep = "")
# axis(side = 1, pos=-.05, lty=0,at = seq(1:6),labels = labs, cex.axis = 1)
#
# dev.off()
#
#
# # D error ~ distance
# pdf(file = "ONE_D_DIST_09.pdf", width = 10, height = 8)
#
# RES=RES9
# boxplot(RES[,c(7:ncol(RES))],varwidth = T,col=c("white","grey"),xaxt="n",
#   ylab="D error rate", xlab="Distance P3/P4",
#   main="D error rate, expected and observed, by distance P3/P4")
# m_res<-apply(RES[,c(7:ncol(RES))],2,function(x) mean(na.omit(x)))
# points(m_res, col='red',pch=18, cex = 2)
# block=1000000/6
# lvl<-levels(cut(1000000, seq(0,1000000,block)))
# vis<-c('ext\n','obs\n')
# vars<-gsub(",","\n",lvl)
# labs<-paste(vis, rep(vars, each = length(vis)), sep = "")
# axis(side = 1,pos=-.05, at = seq(1:length(c(7:ncol(RES)))),labels = labs, cex.axis = 0.75, lty=0)
#
# dev.off()
#
# dir<-list.files(path = ".",pattern = 'samp')
# RES<-as.data.frame(do.call(rbind,lapply(list.files(path = dir[1],
#    pattern = 'test'),function(x) unlist(do_all(paste(dir[1],x,sep="/"))))))
#
# which(RES$Derr_exp == max(RES$Derr_exp))
# plot(RES$Derr_exp)
#
#



################################################################################
################################################################################
################################################################################

# D error ~ sampling

library('parallel')
files<-list.files(path = "sampling/",pattern = 'test')

err_D2<-function(temp){
  temp$FDR<-p.adjust(temp$Pvalue,"fdr")
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
  res<-c("Derr_exp"=D_err_exp_H1,"Derr_obs"=sD_err_exp_H1)
  # res<-c("Derr_exp"=D_err_exp_H1,"DH2err_exp"=D_err_exp_H2,"DH3err_exp"=D_err_exp_H3,"DH123err_exp"=D_err_exp_H123,
  #   "Derr_obs"=sD_err_exp_H1,"DH2err_obs"=sD_err_exp_H2,"DH3err_obs"=sD_err_exp_H3,"DH123err_obs"=sD_err_exp_H123
  # )
  return(res)
}

sampling_error<-function(d){
  # comb<-sample(ind ,d ,replace = FALSE)
  comb<-t(replicate(100,sample(ind,d,replace = FALSE)))
  res<-t(apply(comb,1,function(x) err_D2(subset(data, P1 %in% x & P2 %in% x & P3 %in% x & P4 %in% x))))
  return(res)
}

sampled_D<-function(x){
  data<-read.table(paste("sampling",x,"data.txt",sep="/"),h=T)
  data$EV<-as.factor(paste(data$Donor,data$Recip, sep=""))
  ind<-sort(unique(c(unique(data$P1),unique(data$P2),unique(data$P3),unique(data$P4))))
  lvl<-seq(4,length(ind),4)
  data$FDR<-p.adjust(data$Pvalue,"fdr")
  res<-do.call(cbind,lapply(lvl,function(x) sampling_error(x)))
  m_res<-apply(res,2,function(x) mean(na.omit(x)))
  sd_res<-apply(res,2,function(x) sd(na.omit(x)))
  var_res<-apply(res,2,function(x) var(na.omit(x)))
  RES=c("mean"=m_res,"sd"=sd_res,"var"=var_res)
  return(RES)
}

# lapply(files, function(x) sampled_D(x))

data<-read.table(paste("sampling",files[13],"data.txt",sep="/"),h=T)
data$EV<-as.factor(paste(data$Donor,data$Recip, sep=""))
ind<-sort(unique(c(unique(data$P1),unique(data$P2),unique(data$P3),unique(data$P4))))
lvl<-seq(4,length(ind),4)
# data$FDR<-p.adjust(data$Pvalue,"fdr")
err_D2(data)


res<-do.call(cbind,mclapply(lvl,function(x) sampling_error(x),mc.cores = 20))

pdf(file = "D_SAMPLING_TREE_13.pdf", width = 8, height = 8)
boxplot(res,col=c('white','grey'),xaxt="n",main="Sampling effet on the D statistique error rate \n (expected and observed)",
  ylab="D error rate", xlab="Number of lineages sampled",ylim=c(0,1))
m_res<-apply(res,2,function(x) mean(na.omit(x)))
points(m_res, col='red',pch=18, cex = 1.2)
vis<-c('exp\n','obs\n')
labs<-paste(vis, rep(lvl, each = length(vis)), sep = "")
axis(side = 1,pos=-.04, at = seq(1:ncol(res)),labels = labs, cex.axis = 0.71, lty=0)
dev.off()


# all variance ?
files<-list.files(path = "sampling/",pattern = 'test')


sampled_D<-function(x){
  data<-read.table(paste("sampling",x,"data.txt",sep="/"),h=T)
  data$EV<-as.factor(paste(data$Donor,data$Recip, sep=""))
  ind<-sort(unique(c(unique(data$P1),unique(data$P2),unique(data$P3),unique(data$P4))))
  lvl<-seq(4,length(ind),4)
  # err_D(data)

  res<-do.call(cbind,lapply(lvl,function(x) sampling_error(x)))
  m_res<-apply(res,2,function(x) mean(na.omit(x)))
  sd_res<-apply(res,2,function(x) sd(na.omit(x)))
  var_res<-apply(res,2,function(x) var(na.omit(x)))
  RES=c("mean"=m_res,"sd"=sd_res,"var"=var_res)
  return(RES)
}

x=files[13]
all_samp<-do.call(cbind,mclapply(files, function(x) sampled_D(x),mc.cores = 20))
all_samp<-as.data.frame(t(all_samp))
# boxplot(all_s,xaxt="n",main="Sampling effet on the D statistique error rate \n (expected and observed)",
#   ylab="", xlab="Number of lineages sampled",ylim=c(0,1)))
# m_res<-apply(all_s,2,function(x) mean(na.omit(x)))
# points(m_res, col='red',pch=18, cex = 1.2)
# axis(side = 1,pos=-.04, at = seq(1:ncol(res)),labels = c('meas','sd','var'), cex.axis = 0.71, lty=0)

# sd(error) variation by lvl
boxplot(all_samp[,c(21:40)],xaxt="n",main="D error sd by level of sampling\n(expected and observed)",
  ylab="D error standar deviation", xlab="Number of lineages sampled")
m_res<-apply(all_samp[,c(21:40)],2,function(x) mean(na.omit(x)))
points(m_res, col='red',pch=18, cex = 1.2)
lvl<-seq(4,40,4)
vis<-c('exp\n','obs\n')
labs<-paste(vis, rep(lvl, each = length(vis)), sep = "")
axis(side = 1,pos=-.04, at = seq(1:ncol(all_samp[,c(21:40)])),labels = labs, lty=0)


# var(error) variation by lvl
boxplot(all_samp[,c(41:60)],xaxt="n",main="D error variance by level of sampling\n(expected and observed)",
  ylab="D error variance", xlab="Number of lineages sampled")
m_res<-apply(all_samp[,c(41:60)],2,function(x) mean(na.omit(x)))
points(m_res, col='red',pch=18, cex = 1.2)
lvl<-seq(4,40,4)
vis<-c('exp\n','obs\n')
labs<-paste(vis, rep(lvl, each = length(vis)), sep = "")
axis(side = 1,pos=-.04, at = seq(1:ncol(all_samp[,c(41:60)])),labels = labs, lty=0)


################################################################################
################################################################################
################################################################################

# Error by events, one tree

dir="one_tree"
files<-list.files(path = dir,pattern = 'test')

Devener<-function(x){
  data<-read.table(paste(dir,x,"data.txt",sep="/"),h=T)
  data$EV<-as.factor(paste(data$Donor,data$Recip, sep=""))
  return(data[,c("D","EV")])
}

D3evener<-function(x){
  data<-read.table(paste(dir,x,"data_D3.txt",sep="/"),h=T)
  data$EV<-as.factor(paste(data$Donor,data$Recip, sep=""))
  return(data[,c("D3","EV")])
}

# for D
test<-do.call(rbind,lapply(files, function(x) Devener(x)))

ttest<-test[test$EV %in% c('P1P3','P3P1','P1P4','P4P1','P2P3','P3P2','P2P4','P4P2','N2P1','N2P2','OP1','OP2'),]
table(as.character(ttest$EV))
sum(table(as.character(ttest$EV)))

pdf(file = "D_BY_EV.pdf", width = 8, height = 8)
boxplot(ttest$D~as.character(ttest$EV),width=table(as.character(ttest$EV)),
  main="Dstatistic by transfer events",col='lightgrey',
  xlab="Transfer events", ylab="D statistic",boxwex=T,cex.axis=0.75)
abline(h=0,col="red")
dev.off()

pdf(file = "DALL_BY_EV.pdf")
boxplot(test$D~as.character(test$EV),
  main="D by transfer events",col='lightgrey',
  xlab="Transfer events", ylab="D")
abline(h=0,col="red")
dev.off()

# no outliner
pdf(file = "D_BY_EV_NO_OUT.pdf", width = 8, height = 8)
boxplot(ttest$D~as.character(ttest$EV),width=table(as.character(ttest$EV)),
  main="Dstatistic by transfer events",col='lightgrey',outline=FALSE,
  xlab="Transfer events", ylab="D statistic",boxwex=T)
abline(h=0,col="red")
dev.off()

pdf(file = "DALL_BY_EV_NO_OUT.pdf")
boxplot(test$D~as.character(test$EV),
  main="D by transfer events",col='lightgrey',outline=FALSE,
  xlab="Transfer events", ylab="D")
abline(h=0,col="red")
dev.off()

# For D3
test<-do.call(rbind,lapply(files, function(x) D3evener(x)))

ttest<-test[test$EV %in% c('P1P3','P3P1','P2P3','P3P2','OP1','OP2'),]

pdf(file = "D3_BY_EV.pdf", width = 8, height = 8)
boxplot(ttest$D3~as.character(ttest$EV),width=table(as.character(ttest$EV)),
  main="D3 by transfer events",col='lightgrey',
  xlab="Transfer events", ylab="D3",boxwex=T)
abline(h=0,col="red")
dev.off()

pdf(file = "D3ALL_BY_EV.pdf", width = 10, height = 8)
boxplot(test$D3~as.character(test$EV),
  main="D3 by transfer events",col='lightgrey',
  xlab="Transfer events", ylab="D3")
abline(h=0,col="red")
dev.off()

pdf(file = "D3_BY_EV_NO_OUT.pdf", width = 8, height = 8)
boxplot(ttest$D3~as.character(ttest$EV),width=table(as.character(ttest$EV)),
  main="D3 by transfer events",col='lightgrey',outline=FALSE,
  xlab="Transfer events", ylab="D3",boxwex=T)
abline(h=0,col="red")
dev.off()


pdf(file = "D3ALL_BY_EV_NO_OUT.pdf", width = 10, height = 8)
boxplot(test$D3~as.character(test$EV),
  main="D3 by transfer events",col='lightgrey',outline=FALSE,
  xlab="Transfer events", ylab="D3")
abline(h=0,col="red")
dev.off()

################################################################################
################################################################################
################################################################################
