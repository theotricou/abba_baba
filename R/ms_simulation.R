#!/usr/bin/env Rscript
# Théo

cat("Starting ms simulation \n")
args = commandArgs(trailingOnly=TRUE)

cat("step 1 : initializing \n")

output=args[1]
source(paste(output,"ms_command.R", sep ="/"))
source(paste(output,"sim_parameters", sep ="/"))


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

D_stat <- function(stat_simulation, q){
  p1 = which(da$nu == q[1])
  p2 = which(da$nu == q[2])
  p3 = which(da$nu == q[3])
  p4 = which(da$nu == q[4])
  abba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_abba(x)))
  baba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  ng_p3 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[3]))], "@")[[1]][2])
  ng_p4 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[4]))], "@")[[1]][2])
  dist_p13_p14 <- ng_p4 - ng_p3
  donor = where_is_Daldo(q[1], q[2], q[3], q[4])
  recip = where_is_Raldo(q[1], q[2], q[3], q[4])
  data <- c("P1" = q[1], "P2" = q[2],
  "P3" = q[3], "P4" = q[4],
  "abba" = abba, "baba" = baba, "D" = D,
  "Pvalue" = binom.test(c(abba, baba), p = 0.5)$p.value,
  "d_N2" = dist_p13_p14,
  "Donor" = donor, "Recip" = recip)
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

D3_v2<-function(matrice, t){
  P1<-da[da$nu == t[1],3]
  P2<-da[da$nu == t[2],3]
  P3<-da[da$nu == t[3],3]
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  D3_donor<-position_D_in_tree(t[1],t[2],t[3])
  D3_recip<-position_R_in_tree(t[1],t[2],t[3])
  data <- c(
    "D3" = as.numeric(D3),
    "Donor" = D3_donor,
    "Recip" = D3_recip)
  return(data)
}

################################################################################

### D stat

tree <- read.tree(file.path(output, "spe_tree"))
extant<-grep("l>",tree$tip.label)
tr<-keep.tip(tree, tree$tip.label[extant])
da<-data.frame(na=tree$tip.label[extant], nu=extant, ge=paste("s",1:length(extant),sep=""))

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

cat("step 2 : uniq \n")

uniq <- lapply(rep, function(x){
  if (ncol(x$seg_sites[[1]]) == 1) {
      return(x$seg_sites[[1]][[1]])}})

single_trees<-unlist(lapply(which(uniq != "NULL"), function(x) rep[[x]]$trees[[1]]))

# cat("Outputting trees from single segregating site locus in: ")
# outfile_s <- paste(output, "trees", sep = "/")
# write(single_trees, file=outfile_s)
# cmd <- paste("zip ",outfile_s,".zip ",outfile_s, sep="")
# system(cmd, wait=T)
# cmd <- paste("rm ",outfile_s, sep="")
# system(cmd, wait=T)

cat("step 2 : selecting sites \n")

sites <- do.call("cbind", uniq)
sites <- sites[, colSums(sites != 0) > 1]
topologiesD <- getquatuors(tr)
topologiesD <- apply(topologiesD,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

cat("step 3 : computing D stat \n")

results <- as.data.frame(t(apply(topologiesD, 1, function(x) D_stat(sites, x)))) # bottleneck

cat("step 4 : output D\n")

outfile_d <- paste(output, "data.txt", sep = "/")
write.table(results, outfile_d, sep = "\t", row.names = F, append = F, quote=F)

### D3

topologiesD3 <- gettrios(tr)
topologiesD3 <- apply(topologiesD3,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
single_trees <- lapply(single_trees, function(x) read.tree(text=x))
ord<-sort(row.names(cophenetic(single_trees[[1]])))
blocks<-seq(1,length(single_trees), by=50000)
fromto<-cbind(blocks,c(blocks[-1]-1, length(single_trees)))
list_submat<-lapply(1:length(blocks),function(x) devil_mambojambo(x))
sum_mat<-Reduce('+',list_submat)

cat("step 5 : computing D3\n")

results<-cbind(topologiesD3,as.data.frame(t(apply(topologiesD3, 1, function(x) D3_v2(sum_mat, x)))))

cat("step 6 : output D3\n")

outfile <- paste(output, "data_D3.txt", sep = "/")
write.table(results, outfile, sep = "\t", row.names = F, append = F, quote=F)

cat("End \n")
