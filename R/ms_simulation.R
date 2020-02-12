#!/usr/bin/env Rscript
# Théo

# args = commandArgs(trailingOnly=TRUE)

source("ms_command.R")
library(phangorn)

is_abba <- function(seg_site) {if (seg_site[1] == seg_site[4] & seg_site[2] == seg_site[3] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

is_baba <- function(seg_site) {if (seg_site[1] == seg_site[3] & seg_site[2] == seg_site[4] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

named <- function(P) {return(which(pop == P)[1])}

D_stat <- function(stat_simulation, quatuor){
  abba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_abba(x)))
  baba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  data <- c("P1" = named(quatuor[1]), "P2" = named(quatuor[2]), "P3" = named(quatuor[3]), "P4" = named(quatuor[4]), "abba" = abba, "baba" = baba, "Dstat" = D)
  return(data)
}

is_H_true <- function(test, tree) {
  a_P1 <- c(which(spnd == named(test[[1]])), Ancestors(tree, which(spnd == named(test[[1]])), type = "all"))
  a_P2 <- c(which(spnd == named(test[[2]])), Ancestors(tree, which(spnd == named(test[[2]])), type = "all"))
  a_P3 <- c(which(spnd == named(test[[3]])), Ancestors(tree, which(spnd == named(test[[3]])), type = "all"))
  first_son_P3 <- a_P3[match(mrca.phylo(tree, c(which(spnd == named(test[[3]])), which(spnd == named(test[[2]])))), a_P3) - 1]
  if (test[[7]] < 0) {
    first_son_P1 <- a_P1[match(mrca.phylo(tree, c(which(spnd == named(test[[1]])), which(spnd == named(test[[2]])))), a_P1) - 1]
    if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
        true_recip %in% a_P1 ||
        true_donor %in% Descendants(tree, first_son_P1, type = "all") &
        true_recip %in% a_P3) {return("True")} else {return("False")}
  } else if (test[[7]] > 0) {
    first_son_P2 <- a_P2[match(mrca.phylo(tree, c(which(spnd == named(test[[1]])), which(spnd == named(test[[2]])))), a_P2) - 1]
    if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
        true_recip %in% a_P2 ||
        true_donor %in% Descendants(tree, first_son_P2, type = "all") &
        true_recip %in% a_P3) {return("True")} else {return("False")}
  } else {return("False")}
}

for (i in 2:length(pop)) {pop[i] <- pop[i-1] + pop[i]}

rep = simulate(model, nsim = 50000, core = 4)
uniq <- lapply(rep, function(x){
  if (ncol(x$seg_sites[[1]][[1]]) == 1) {
      return(x$seg_sites[[1]][[1]])}})

sites <- c()
for (i in 1:length(uniq)) {
  if (length(uniq[[i]]) != 0) {
    sites <- cbind(sites, uniq[[i]])}}

tree <- read.tree("Species_tree")
spnd<-c(tree$tip.label, tree$node.label)
dist_to_root <- dist.nodes(tree)[,7]

l1=l2=l3=l4 = seq(from = 1, to = nrow(sites))
library("data.table")
comb <- t(CJ(l1,l2,l3,l4)[l1 != l2 & l1 != l3 & l1 != l4 & l2 != l3 & l2 != l4 & l3 != l4])
rem <- c()
for (i in 1:ncol(comb)) {
  if (! dist_to_root[mrca.phylo(tree, c(named(comb[1,i]), named(comb[2,i])))] > dist_to_root[mrca.phylo(tree, c(named(comb[2,i]), named(comb[3,i])))] ||
  ! dist_to_root[mrca.phylo(tree, c(named(comb[2,i]), named(comb[3,i])))] > dist_to_root[mrca.phylo(tree, c(named(comb[3,i]), named(comb[4,i])))]) {
    rem <- c(rem, i)
}}
comb_f = t(unique(apply(comb[,-rem], 2, function(x) {return(c(min(x[c(1,2)]), max(x[c(1,2)]), x[3], x[4]))})))
topologies <- t(unique(comb_f[duplicated(comb_f),]))
results <- as.data.frame(t(apply(topologies, 2, function(x) D_stat(sites, x))))

# now test if basic H1 is true: P3 introgressed in P2 (D > 0 or P1 (D < 1)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file
results$IsH1 <- apply(results, 1, function(x) is_H_true(x ,tree))
print(results)




#
#
#
# for (i in 1:nrow(results)){
# test = results[i,]
#   a_P1 <- c(which(spnd == named(test[[1]])), Ancestors(tree, which(spnd == named(test[[1]])), type = "all"))
#   a_P2 <- c(which(spnd == named(test[[2]])), Ancestors(tree, which(spnd == named(test[[2]])), type = "all"))
#   a_P3 <- c(which(spnd == named(test[[3]])), Ancestors(tree, which(spnd == named(test[[3]])), type = "all"))
#   first_son_P3 <- a_P3[match(mrca.phylo(tree, c(which(spnd == named(test[[3]])), which(spnd == named(test[[2]])))), a_P3) - 1]
#   if (test[[7]] < 0) {
#     first_son_P1 <- a_P1[match(mrca.phylo(tree, c(which(spnd == named(test[[1]])), which(spnd == named(test[[2]])))), a_P1) - 1]
#     if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
#         true_recip %in% a_P1 ||
#         true_donor %in% Descendants(tree, first_son_P1, type = "all") &
#         true_recip %in% a_P3) {print("True")} else {print("False")}
#   } else if (test[[7]] > 0) {
#     first_son_P2 <- a_P2[match(mrca.phylo(tree, c(which(spnd == named(test[[1]])), which(spnd == named(test[[2]])))), a_P2) - 1]
#     if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
#         true_recip %in% a_P2 ||
#         true_donor %in% Descendants(tree, first_son_P2, type = "all") &
#         true_recip %in% a_P3) {print("True")} else {print("False")}
#   } else {print("False")}
# }
#
# apply(results, 1, function(x) is_H_true(x ,tree))
