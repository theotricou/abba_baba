#!/usr/bin/env Rscript
# Théo

args = commandArgs(trailingOnly=TRUE)
output=args[1]

output="test669"
source(paste(output,"ms_command.R", sep ="/"))

spe<-read.tree(paste(output,"spe_tree", sep = "/"))
extant<-grep("@0",spe$tip.label)
names<-data.frame(a=paste("s",1:length(extant),sep=""), b=paste(grep("@0",spe$tip.label), "@0",sep=""))
extant<-keep.tip(spe, spe$tip.label[extant])

validateandreorder<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  return(names(sort(apply(submat,1,sum))))
}

gettrios<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  alltrios<-combn(tr$tip.label,3)
  RES<-t(apply(alltrios, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}

topologies <- gettrios(extant)

trees_file<-unz(paste(output,"trees.zip", sep ="/"),paste(output,"trees", sep ="/"))
trees = read.tree(trees_file)
close(trees_file)

devil_mambojambo<-function(range){
return(Reduce('+', lapply(trees[as.integer(fromto[range,1]):as.integer(fromto[range,2])], function(x) cophenetic(x)[ord,ord])))
}

nodes<-c(spe$tip.label, spe$node.label)
donor<-which(nodes == name_donor)
recip<-which(nodes == name_recip)
ord=sort(row.names(cophenetic(trees[[1]])))

blocks<-seq(1,length(trees), by=1000)
fromto<-cbind(blocks,c(blocks[-1]-1, length(trees)))
list_submat<-lapply(1:length(blocks),function(x) devil_mambojambo(x))
sum_mat<-Reduce('+',list_submat)

leaf_descendants <- function(leaf, comp_leaf) {
  all_anc_leaf <- c(leaf, Ancestors(spe, leaf, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(spe, c(leaf, comp_leaf)), all_anc_leaf) - 1]
  des_leaf <- Descendants(spe, as.integer(first_son_leaf), type = "all")
  return(des_leaf)
}

leaf_ancstors<-function(p, p_ext){
  all_anc <- c(p, Ancestors(spe, p, type = "all"))
  anc_only <- all_anc[1:match(mrca.phylo(spe, c(p, p_ext)), all_anc) - 1]
  return(anc_only)
}

position_D_in_tree<-function(p1,p2,p3){
  P3des<-leaf_descendants(p3,p2)
  P2des<-leaf_descendants(p2,p1)
  P1des<-leaf_descendants(p1,p2)
  if (donor %in% P3des){return("P3")}
  else if (donor %in% P2des){return("P2")}
  else if (donor %in% P1des){return("P1")}
  else if (donor %in% Descendants(spe,mrca.phylo(spe, c(p1, p3)), type="all")){return("N1")}
  else{return("O")}
}

position_R_in_tree<-function(p1,p2,p3){
  P3anc<-leaf_ancstors(p3,p2)
  P2anc<-leaf_ancstors(p2,p1)
  P1anc<-leaf_ancstors(p1,p2)
  if (recip %in% P3anc){return("P3")}
  else if (recip %in% P2anc){return("P2")}
  else if (recip %in% P1anc){return("P1")}
  else if (recip %in% Descendants(spe,mrca.phylo(spe, c(p1, p3)), type="all")){return("N1")}
  else{return("O")}
}
D3_v2<-function(matrice, trio){
  P1<-as.character(names[names$b == trio[1],1])
  P2<-as.character(names[names$b == trio[2],1])
  P3<-as.character(names[names$b == trio[3],1])
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  p1=which(nodes == trio[1])
  p2=which(nodes == trio[2])
  p3=which(nodes == trio[3])
  D3_donor<-position_D_in_tree(p1,p2,p3)
  D3_recip<-position_R_in_tree(p1,p2,p3)
  data <- c(
    "D3" = as.numeric(D3),
    "Donor" = D3_donor,
    "Recip" = D3_recip)
  return(data)
}

results<-cbind(topologies,as.data.frame(t(apply(topologies, 1, function(x) D3_v2(sum_mat, x)))))
outfile <- paste(output, "data_D3.txt", sep = "/")
write.table(results, outfile, sep = "\t", row.names = F, append = F, quote=F)
