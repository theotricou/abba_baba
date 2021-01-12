#!/usr/bin/env Rscript
# Theo
# DeeD_reader

library("phangorn")
library('reshape2')

get_neighbour<-function(p,r){
  all_anc<-c(p,Ancestors(t,p,type = "all"))
  anc<-c(all_anc[1:match(mrca.phylo(t, c(p, r)), all_anc) - 1])
  last_anc<-tail(anc, n=1)
  desc_P<-c(last_anc,Descendants(t,last_anc, type = "all"))
  res<-list("anc"=unique(anc),
    "sis"=unique(desc_P)
    )
  return(res)
}

get_neighbour_and_more<-function(p,r,q){
  # this gets you N2 branches
  p123<-mrca.phylo(t, c(r,q))
  p1234<-mrca.phylo(t, c(p,r))
  all_anc_n2<-c(p123,Ancestors(t,p123, type='all'))
  n2<-all_anc_n2[1:match(p1234,all_anc_n2) - 1]
  bis<-c(tail(n2, n=1),Descendants(t,tail(n2, n=1),type = "all"))
  p123_d=c(Descendants(t,p123, type="all"))
  n2=setdiff(bis,p123_d)
  # this gets you all P4 and outsider
  all<-seq(1,length(nodes))
  desc_P<-setdiff(all,bis)
  # and this all P4 ancestor
  all_anc<-c(p,Ancestors(t,p,type = "all"))
  anc<-c(all_anc[1:match(p1234, all_anc) - 1])
  res<-list("anc"=unique(anc),
    "sis"=unique(desc_P),
    "n2"=unique(n2)
    )
  return(res)
}

P_sistor_ancestor<-function(quat){
  p1<-get_neighbour(quat[1],quat[2])
  p2<-get_neighbour(quat[2],quat[1])
  p3<-get_neighbour(quat[3],quat[1])
  p4<-get_neighbour_and_more(quat[4],quat[1],quat[3])
  res<-list(
    "sP1"<-p1$sis,
    "aP1"<-p1$anc,
    "sP2"<-p2$sis,
    "aP2"<-p2$anc,
    "sP3"<-p3$sis,
    "aP3"<-p3$anc,
    "sP4"<-p4$sis,
    "aP4"<-p4$anc,
    "N2"<-p4$n2
  )
  return(res)
}

plusor<-function(x){
  # more ABBA
  quat=c(as.numeric(as.character(x[1:4])))
  block=P_sistor_ancestor(quat)
  # P1->P4
  a<-expand.grid(block[[1]],block[[8]])
  # P4->P1
  b<-expand.grid(block[[7]],block[[2]])
  # P2->P3
  c<-expand.grid(block[[3]],block[[6]])
  # P3->P2
  d<-expand.grid(block[[5]],block[[4]])
  # N2->P1
  e<-expand.grid(block[[9]],block[[2]])
  res<-rbind(a,b,c,d,e)
  res2<-apply(res, c(1,2), function(x) nodes[x])
  plus<-apply(res2, 1, function(x) mat[x[1],x[2]]<<-mat[x[1],x[2]]+1)
  return(res2)
}

minusor<-function(x){
  # more BABA
  quat=c(as.numeric(as.character(x[1:4])))
  block=P_sistor_ancestor(quat)
  # P1->P3
  a<-expand.grid(block[[1]],block[[6]])
  # P3->P1
  b<-expand.grid(block[[5]],block[[2]])
  # P2->P4
  c<-expand.grid(block[[3]],block[[8]])
  # P4->P2
  d<-expand.grid(block[[7]],block[[4]])
  # N2->P2
  e<-expand.grid(block[[9]],block[[4]])
  res<-rbind(a,b,c,d,e)
  res2<-apply(res, c(1,2), function(x) nodes[x])
  plus<-apply(res2, 1, function(x) mat[x[1],x[2]]<<-mat[x[1],x[2]]+1)
  return(res2)
}

nullor<-function(x){
  # no excess
  quat=c(as.numeric(as.character(x[1:4])))
  block=P_sistor_ancestor(quat)
  # P1->P3
  a1<-expand.grid(block[[1]],block[[6]])
  # P1->P4
  a2<-expand.grid(block[[1]],block[[8]])
  # P2->P3
  c1<-expand.grid(block[[3]],block[[6]])
  # P2->P4
  c2<-expand.grid(block[[3]],block[[8]])
  # P3->P1
  b1<-expand.grid(block[[5]],block[[2]])
  # P3->P2
  d1<-expand.grid(block[[5]],block[[4]])
  # P4->P1
  b2<-expand.grid(block[[7]],block[[2]])
  # P4->P2
  d2<-expand.grid(block[[7]],block[[4]])
  # N2->P2
  e1<-expand.grid(block[[9]],block[[4]])
  # N2->P1
  e2<-expand.grid(block[[9]],block[[2]])
  res<-rbind(a1,b1,c1,d1,e1,a2,b2,c2,d2,e2)
  res2<-apply(res, c(1,2), function(x) nodes[x])
  minus<-apply(res2, 1, function(x) mat[x[1],x[2]]<<-mat[x[1],x[2]]-1)
  return(res2)
}

eventer<-function(x){
  if(as.numeric(x)[1]<0.05){
    if(as.numeric(x[2])>0){
      return("+")
    }else{
      return("-")
    }
  }else{
    return("0")
  }
}

D_reader<-function(x){
  if(x[19]=="+"){
    plusor(x)
  }else if(x[19]=="-"){
    minusor(x)
  }else
  if(x[19]=="0"){
    nullor(x)
  }
}

d<-read.table('data.txt', h=T) # format P1, P2, P3, P4, abba, baba, D, Pvalue
tree<-read.tree("spe_tree")
extant_nodes<-c(tree$tip.label, tree$node.label)
extant<-grep("ali>",tree$tip.label)
t<-keep.tip(tree, tree$tip.label[extant])


nodes<-c(t$tip.label, t$node.label)
dict_nodes<-as.data.frame(cbind(nodes, seq(1,length(nodes))))
dict_nodes$V3<-lapply(nodes, function(x) which(x == extant_nodes))

d[,c(1,2,3,4)] <- apply(d[,c(1,2,3,4)], c(1,2), function(x) return(dict_nodes[which(dict_nodes$V3 == x),2]))

d$fdr<-as.numeric(p.adjust(d$Pvalue,"fdr"))
d$event<-as.factor(apply(d,1,function(x) eventer(x[c(18,7)])))

# create matrix of event, row = donor, col = recipient.
mat<-matrix(0, nrow = length(nodes), ncol = length(nodes))
colnames(mat)<-nodes
rownames(mat)<-nodes
res<-apply(d, 1, function(x) D_reader(x))
bb=melt(mat)
print(bb[order(bb$value),])

#
#
# check<-function(x){
#   if (class(x) == "matrix"){
#     for (i in 1:nrow(x)){
#       if (x[i,1] %in% c("ali>8@0", "ali>7@0")){
#         if (x[i,2] %in% c("ali>8@0", "ali>7@0")){
#           return(x)
#         }
#       }
#     }
#   }
# }
#
# lapply(res, function(x) check(x))
#
#
#
# true_eventer<-function(x){
#   if(paste(x[1],x[2],sep="")%in%c("P2P3","P1P4","P3P2","P4P1","N2P1")){
#     return("+")
#   }else if (paste(x[1],x[2],sep="")%in%c("P1P3","P2P4","P3P1","P4P2","N2P2")){
#     return("-")
#   }else{
#     return("0")
#   }
# }
#
# d$tevent<-as.factor(apply(d,1,function(x) true_eventer(x[c(9,10)])))
# D_true_reader<-function(x){
#   if(x[20]=="+"){
#     plusor(x)
#   }else if(x[20]=="-"){
#     minusor(x)
#   }else{
#     nullor(x)
#   }
# }
#
# mat<-matrix(0, nrow = length(nodes), ncol = length(nodes))
# colnames(mat)<-nodes
# rownames(mat)<-nodes
# res<-apply(d, 1, function(x) D_true_reader(x))
# mat
#
# bb=melt(mat)
# bb[order(bb$value),]
#
