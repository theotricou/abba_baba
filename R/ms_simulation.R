#!/usr/bin/env Rscript
# Théo

cat("Starting ms simulation \n")
args = commandArgs(trailingOnly=TRUE)

cat("step 1 : initializing \n")

output=args[1]
# output="."
source(paste(output,"ms_command.R", sep ="/"))
source(paste(output,"sim_parameters", sep ="/"))
library(parallel)


is_abba <- function(seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[4]){return(1)} else {return(0)}}
is_baba <- function(seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[3]){return(1)} else {return(0)}}

validateandreorderD<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
}}

getquatuors<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorderD(x, dist)))
  return(RES)
}

leaf_descendants <- function(p, p_ext) {
  all_anc_leaf <- c(p, Ancestors(tree, p, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(tree, c(p, p_ext)), all_anc_leaf) - 1]
  des_leaf <- c(first_son_leaf,Descendants(tree, as.integer(first_son_leaf), type = "all"))
  return(des_leaf)
}

leaf_ancestors<-function(p, p_ext){
  all_anc <- c(p, Ancestors(tree, p, type = "all"))
  anc_only <- all_anc[1:match(mrca.phylo(tree, c(p, p_ext)), all_anc) - 1]
  return(anc_only)
}

where_is_Raldo<-function(p1,p2,p3,p4){
  if (true_recip %in% leaf_ancestors(p1,p2)) { return("P1")}
  else if (true_recip %in% leaf_ancestors(p2,p1)) { return("P2")}
  else if (true_recip %in% leaf_ancestors(p3,p1)) { return("P3")}
  else if (true_recip %in% leaf_ancestors(p4,p1)) { return("P4")}
  else {
    if (!true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")) {return("O")}
    else {
      if (true_recip %in% leaf_ancestors(p1,p3)) {return("N1")}
      else {
        if (true_recip %in% leaf_ancestors(p1,p4)) {return("N2")}
        else {return("O")}
}}}}

where_is_Daldo<-function(p1,p2,p3,p4){
  if (true_donor %in% leaf_descendants(p1,p2)) { return("P1")}
  else if (true_donor %in% leaf_descendants(p2,p1)) { return("P2")}
  else if (true_donor %in% leaf_descendants(p3,p1)) { return("P3")}
  else if (true_donor %in% leaf_descendants(p4,p1)) { return("P4")}
  else {
    if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")) {return("O")}
    else {
      if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")) {return("N2")}
      else {return("N1")}
}}}

N2_diversity<-function(q){
  a<-mrca.phylo(tree, spnd[unlist(q)])
  p1p4Desc<-Descendants(tree, mrca.phylo(tree, spnd[q]),type="tip")[[1]]
  p1p3Desc<-Descendants(tree, mrca.phylo(tree, spnd[q][c(1:3)]),type="tip")[[1]]
  p0<-c(q[4],Ancestors(tree, q[4], type = "all"))
  p4Desc<-Descendants(tree,p0[which(p0 == a)-1],type="tip")[[1]]
  all<-c(p1p3Desc,p4Desc)
  `%notin%` <- Negate(`%in%`)
  N2_tips<-p1p4Desc[p1p4Desc %notin% all]
  n_hidden<-length(N2_tips)
  n_ext<-length(grep("s>",tree$tip.label[N2_tips]))
  return(c("H"=n_hidden, "E"=n_ext))
}

D_stat <- function(stat_simulation, q){
  p1 = which(da$nu == q[1])
  p2 = which(da$nu == q[2])
  p3 = which(da$nu == q[3])
  p4 = which(da$nu == q[4])
  abba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_abba(x)))
  baba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  donor = where_is_Daldo(q[1], q[2], q[3], q[4])
  recip = where_is_Raldo(q[1], q[2], q[3], q[4])
  N2_tips<-N2_diversity(c(q[1], q[2], q[3], q[4]))
  data <- c(
    "P1" = q[1],
    "P2" = q[2],
    "P3" = q[3],
    "P4" = q[4],
    "abba" = abba,
    "baba" = baba,
    "D" = D,
    "Pvalue" = binom.test(c(abba, baba), p = 0.5)$p.value,
    "Donor" = donor,
    "Recip" = recip,
    "dP1P2" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[2]))], "@")[[1]][2]),
    "dP1P3" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[3]))], "@")[[1]][2]),
    "dP1P4" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[4]))], "@")[[1]][2]),
    "N2_H"=N2_tips[1],
    "N2_E"=N2_tips[2]
  )
  return(data)
}

D3_v2<-function(matrice, t){
  P1<-da[da$nu == t[1],3]
  P2<-da[da$nu == t[2],3]
  P3<-da[da$nu == t[3],3]
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  D3_donor<-position_D_in_tree(t[1],t[2],t[3])
  D3_recip<-position_R_in_tree(t[1],t[2],t[3])
  data <- c(
    "P1" = t[1],
    "P2" = t[2],
    "P3" = t[3],
    "dp23" = matrice[P2,P3],
    "dp13" = matrice[P1,P3],
    "D3" = as.numeric(D3),
    "Donor" = D3_donor,
    "Recip" = D3_recip,
    "dP1P2" = as.integer(strsplit(spnd[mrca.phylo(tree, c(t[1], t[2]))], "@")[[1]][2]),
    "dP1P3" = as.integer(strsplit(spnd[mrca.phylo(tree, c(t[1], t[3]))], "@")[[1]][2])
  )
  return(data)
}

validateandreorderD3<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  return(names(sort(apply(submat,1,sum))))
}

gettrios<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  alltrios<-combn(tr$tip.label,3)
  RES<-t(apply(alltrios, 2, function(x) validateandreorderD3(x, dist)))
  return(RES)
}

devil_mambojambo<-function(range){
  return(Reduce('+', lapply(single_trees[as.integer(fromto[range,1]):as.integer(fromto[range,2])], function(x) cophenetic(x)[ord,ord])))
}

leaf_descendants <- function(leaf, comp_leaf) {
  all_anc_leaf <- c(leaf, Ancestors(tree, leaf, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(tree, c(leaf, comp_leaf)), all_anc_leaf) - 1]
  des_leaf <- c(first_son_leaf,Descendants(tree, as.integer(first_son_leaf), type = "all"))
  return(des_leaf)
}

leaf_ancestors<-function(p, p_ext){
  all_anc <- c(p, Ancestors(tree, p, type = "all"))
  anc_only <- all_anc[1:match(mrca.phylo(tree, c(p, p_ext)), all_anc) - 1]
  return(anc_only)
}

position_D_in_tree<-function(p1,p2,p3){
  if (true_donor %in% leaf_descendants(p3,p2)){return("P3")}
  else if (true_donor %in% leaf_descendants(p2,p1)){return("P2")}
  else if (true_donor %in% leaf_descendants(p1,p2)){return("P1")}
  else if (true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){return("N1")}
  else{return("O")}
}

position_R_in_tree<-function(p1,p2,p3){
  if (true_recip %in% leaf_ancestors(p3,p2)){return("P3")}
  else if (true_recip %in% leaf_ancestors(p2,p1)){return("P2")}
  else if (true_recip %in% leaf_ancestors(p1,p2)){return("P1")}
  else if (true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){return("N1")}
  else{return("O")}
}

uniquer<-function(x){
  if (ncol(x$seg_sites[[1]]) == 1) { # ajouter filtre sur site sum(seg) == 1
      return(x$seg_sites[[1]][[1]])}
}

is_sister<-function(x){
  donor<-as.character(unlist(x[9]))
  if (donor %in% c("O",'N2','N1')){
    return('Nope')
  }else{
    col<-as.numeric(strsplit(donor,'P')[[1]][2])
    node=as.numeric(unlist(x[col]))
    if (true_donor %in% c(node, Ancestors(tree, node))){
      return('True')
    }else{
      return('Sister')
    }
  }
}

################################################################################

### D stat

tree <- read.tree(file.path(output, "spe_tree"))
extant<-grep("l>",tree$tip.label)
tr<-keep.tip(tree, tree$tip.label[extant])
da<-data.frame(na=tree$tip.label[extant], nu=extant, ge=paste("s",1:length(extant),sep=""))
# max_d_to_extant_outgroup <- max(cophenetic(tr))/2

cat("step 2 : running simulation \n")

if (SEED == 0) {rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
} else {
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  M <- N_CORE
  s <- .Random.seed
  for (i in 1:M){s <- nextRNGStream(s)}
  rep <- mclapply(X = 1:N_SIMULATION,
               FUN = function(x) simulate(model),
               mc.cores = M,
               mc.set.seed = TRUE
  )
}

uniq<-mclapply(rep, function(x) uniquer(x),mc.cores = N_CORE)

single_trees<-unlist(lapply(which(uniq != "NULL"), function(x) rep[[x]]$trees[[1]]))

# cat("Outputting trees from single segregating site locus in: ")
# outfile_s <- paste(output, "trees", sep = "/")
# write(single_trees, file=outfile_s)
# cmd <- paste("zip ",outfile_s,".zip ",outfile_s, sep="")
# system(cmd, wait=T)
# cmd <- paste("rm ",outfile_s, sep="")
# system(cmd, wait=T)

sites <- do.call("cbind", uniq)
sites <- sites[, colSums(sites != 0) > 1]
topologiesD <- getquatuors(tr)
topologiesD <- apply(topologiesD,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

cat("step 3 : computing D stat \n")

results<-as.data.frame(do.call(rbind,mclapply(as.data.frame(t(topologiesD)), function(x) D_stat(sites, x),mc.cores = N_CORE)))

results$is_sister<-apply(results,1, function(x) is_sister(x))

cat("step 4 : output D\n")

outfile_d <- paste(output, "data.txt", sep = "/")
write.table(results, outfile_d, sep = "\t", row.names = F, append = F, quote=F)

### D3
#
# topologiesD3 <- gettrios(tr)
# topologiesD3 <- apply(topologiesD3,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
# colnames(topologiesD3)<-c('P1','P2','P3')
# single_trees <- lapply(single_trees, function(x) read.tree(text=x))
# ord<-sort(row.names(cophenetic(single_trees[[1]])))
# blocks<-seq(1,length(single_trees), by=50000)
# fromto<-cbind(blocks,c(blocks[-1]-1, length(single_trees)))
# list_submat<-lapply(1:length(blocks),function(x) devil_mambojambo(x))
# sum_mat<-Reduce('+',list_submat)
#
# cat("step 5 : computing D3\n")
#
# results<-do.call(rbind,mclapply(as.data.frame(t(topologiesD3)), function(x) D3_v2(sum_mat, x),mc.cores = N_CORE))

cat("step 6 : output D3\n")

outfile <- paste(output, "data_D3.txt", sep = "/")
write.table(results, outfile, sep = "\t", row.names = F, append = F, quote=F)

cat("End \n")
