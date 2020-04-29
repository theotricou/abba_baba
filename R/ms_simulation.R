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
      }
    }
  }
}

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
    }
  }
}

D_stat <- function(stat_simulation, quatuor){
  P1 = as.integer(strsplit(quatuor[1], "@")[[1]][1])
  p1 = pop[P1]
  P2 = as.integer(strsplit(quatuor[2], "@")[[1]][1])
  p2 = pop[P2]
  P3 = as.integer(strsplit(quatuor[3], "@")[[1]][1])
  p3 = pop[P3]
  P4 = as.integer(strsplit(quatuor[4], "@")[[1]][1])
  p4 = pop[P4]
  abba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_abba(x)))
  baba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  ng_p3 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(P1, P3))], "@")[[1]][2])
  ng_p4 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(P1, P4))], "@")[[1]][2])
  dist_p13_p14 <- ng_p4 - ng_p3
  donor = where_is_Daldo(P1, P2, P3, P4)
  recip = where_is_Raldo(P1, P2, P3, P4)
  data <- c("P1" = quatuor[1], "P2" = quatuor[2],
  "P3" = quatuor[3], "P4" = quatuor[4],
  "abba" = abba, "baba" = baba, "D" = D,
  "Pvalue" = binom.test(c(abba, baba), p = 0.5, conf.level= 0.99)$p.value,
  "d_N2" = dist_p13_p14,
  "Donor" = donor, "Recip" = recip)
  return(data)
}

tree <- read.tree(file.path(output, "spe_tree"))
tr = keep.tip(tree, tree$tip.label[which(pop == 1)])

for (i in 2:length(pop)) {pop[i] <- pop[i-1] + pop[i]}

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

cat("Outputting trees from single segregating site locus in: ")
single_trees=unlist(lapply(which(uniq != "NULL"), function(x) rep[[x]]$trees[[1]]))
outfile_s <- paste(output, "trees", sep = "/")
write(single_trees, file=outfile_s)
cmd <- paste("zip ",outfile_s,".zip ",outfile_s, sep="")
system(cmd, wait=T)
cmd <- paste("rm ",outfile_s, sep="")
system(cmd, wait=T)

cat("step 2 : sites \n")

sites <- do.call("cbind", uniq)
sites <- sites[, colSums(sites != 0) > 1]

topologies <- getquatuors(tr)
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

cat("step 3 : computing summary statistics \n")

results <- as.data.frame(t(apply(topologies, 1, function(x) D_stat(sites, x)))) # bottleneck

cat("step 4 : output D\n")

outfile_d <- paste(output, "data.txt", sep = "/")
write.table(results, outfile_d, sep = "\t", row.names = F, append = F, quote=F)
