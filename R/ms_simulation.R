#!/usr/bin/env Rscript
# Théo

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  if(file.exists("Simulation/ms_command.R")) {
    source("Simulation/ms_command.R")
    output = "Simulation"
  } else {
    print("You need to creat a ms command by running build_ms_command.py to run this script")
    quit()
  }
} else if (length(args) == 1) {
  source(args[1])
  path = strsplit(args[1], "/")[[1]]
  if (length(path) == 1) {output = "./"} else {output = paste(path[seq(1,length(path) -1)], collapse = '/')}
} else if (length(args) == 2) {
  source(args[1])
  output = args[2]
} else {
  print("Too much argument were given, only the two first one will be used as input and output respectively")
  source(args[1])
  path = strsplit(args[2], "/")[[1]]
  if (length(path) == 1) {output = "./"} else {output = paste(path[seq(1,length(path) -1)], collapse = '/')}
}

source("Simulation/ms_command.R")
output = "Simulation"

is_abba <- function(seg_site) {if (seg_site[1] == seg_site[4] & seg_site[2] == seg_site[3] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

is_baba <- function(seg_site) {if (seg_site[1] == seg_site[3] & seg_site[2] == seg_site[4] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

named <- function(P) {return(which(pop == P)[1])}

D_stat <- function(stat_simulation, quatuor){
  abba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_abba(x)))
  baba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  data <- c("P1" = named(quatuor[1]), "P2" = named(quatuor[2]), "P3" = named(quatuor[3]), "P4" = named(quatuor[4]),
   "abba" = abba, "baba" = baba, "Dstat" = D, "Pvalue" = binom.test(c(abba, baba), p = 0.5, conf.level= 0.99)$p.value)
  return(data)
}

for (i in 2:length(pop)) {pop[i] <- pop[i-1] + pop[i]}

rep = simulate(model, nsim = 500000, cores = 3, seed = 23805)
# rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE, seed = SEED)
uniq <- lapply(rep, function(x){
  if (ncol(x$seg_sites[[1]][[1]]) == 1) {
      return(x$seg_sites[[1]][[1]])}})

sites <- c()
for (i in 1:length(uniq)) {
  if (length(uniq[[i]]) != 0 & sum(uniq[[i]]) > 1) {
    sites <- cbind(sites, uniq[[i]])}}

tree <- read.tree(file.path(output, "spe_tree"))
spnd<-c(tree$tip.label, tree$node.label)
dist_to_root <- dist.nodes(tree)[,length(tree$tip.label)+1]

l1=l2=l3=l4 = seq(from = 1, to = nrow(sites))
comb <- t(CJ(l1,l2,l3,l4)[l1 != l2 & l1 != l3 & l1 != l4 & l2 != l3 & l2 != l4 & l3 != l4])
rem <- c()
for (i in 1:ncol(comb)) {
  p1 <- named(comb[1,i])
  p2 <- named(comb[2,i])
  p3 <- named(comb[3,i])
  p4 <- named(comb[4,i])
  d12 <- dist_to_root[mrca.phylo(tree, c(p1, p2))]
  d13 <- dist_to_root[mrca.phylo(tree, c(p1, p3))]
  d14 <- dist_to_root[mrca.phylo(tree, c(p1, p4))]
  if (d12 > d13 & d12 > d14 & d13 > d14) {
    rem <- c(rem, i)
}}
comb_f = t(unique(apply(comb[,rem], 2, function(x) {return(c(min(x[c(1,2)]), max(x[c(1,2)]), x[3], x[4]))})))
topologies <- t(unique(comb_f[duplicated(comb_f),]))
results <- as.data.frame(t(apply(topologies, 2, function(x) D_stat(sites, x))))

true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file
# results$IsH1 <- apply(results, 1, function(x) is_H_true(x ,tree))


# level 0 = bug, level 1 no signal, level 2 no signal but should, level 3 H1 okayed, level 4 H2 okayed
results$H0 <- 0
for (i in 1:nrow(results)){
test = results[i,]
a_P1 <- c(test[[1]], Ancestors(tree, test[[1]], type = "all"))
a_P2 <- c(test[[2]], Ancestors(tree, test[[2]], type = "all"))
a_P3 <- c(test[[3]], Ancestors(tree, test[[3]], type = "all"))
a_P4 <- c(test[[4]], Ancestors(tree, test[[4]], type = "all"))
first_son_P1 <- a_P1[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P1) - 1]
first_son_P2 <- a_P2[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P2) - 1]
first_son_P3 <- a_P3[match(mrca.phylo(tree, c(test[[1]], test[[3]])), a_P3) - 1]
first_son_P4 <- a_P4[match(mrca.phylo(tree, c(test[[1]], test[[4]])), a_P4) - 1]

if (test[[7]] < 0 & test[[8]] < 0.01) {
  if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
      true_recip %in% a_P1 ||
      true_donor %in% Descendants(tree, first_son_P1, type = "all") &
      true_recip %in% a_P3) {results[i,9] <- 3} else {results[i,9] <- 4}
} else if (test[[7]] > 0 & test[[8]] < 0.01) {
  if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
      true_recip %in% a_P2 ||
      true_donor %in% Descendants(tree, first_son_P2, type = "all") &
      true_recip %in% a_P3) {results[i,9] <- 3} else {results[i,9] <- 4}
} else if (test[[8]] >= 0.01) {
  if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
      true_recip %in% a_P1 ||
      true_donor %in% Descendants(tree, first_son_P1, type = "all") &
      true_recip %in% a_P3 ||
      true_donor %in% Descendants(tree, first_son_P3, type = "all") &
      true_recip %in% a_P2 ||
      true_donor %in% Descendants(tree, first_son_P2, type = "all") &
      true_recip %in% a_P3) {results[i,9] <- 2} else {results[i,9] <- 1}
  }
}

for (i in 1:nrow(results)){
  test = results[i,]
  a_P1 <- c(test[[1]], Ancestors(tree, test[[1]], type = "all"))
  a_P2 <- c(test[[2]], Ancestors(tree, test[[2]], type = "all"))
  a_P3 <- c(test[[3]], Ancestors(tree, test[[3]], type = "all"))
  a_P4 <- c(test[[4]], Ancestors(tree, test[[4]], type = "all"))
  first_son_P1 <- a_P1[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P1) - 1]
  first_son_P2 <- a_P2[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P2) - 1]
  first_son_P3 <- a_P3[match(mrca.phylo(tree, c(test[[1]], test[[3]])), a_P3) - 1]
  first_son_P4 <- a_P4[match(mrca.phylo(tree, c(test[[1]], test[[4]])), a_P4) - 1]
  d_P1 <- Descendants(tree, first_son_P1, type = "all")
  d_P2 <- Descendants(tree, first_son_P2, type = "all")
  d_P3 <- Descendants(tree, first_son_P3, type = "all")
  d_P4 <- Descendants(tree, first_son_P4, type = "all")
  d_quatuor = Descendants(tree, mrca.phylo(tree, c(test[[1]], test[[4]])), type = "all")
  if (test[[8]] >= 0.05) {
    if (true_donor %in% d_P1 & true_recip %in% a_P3 ||
    true_donor %in% d_P2 & true_recip %in% a_P3 ||
    true_donor %in% d_P3 & true_recip %in% a_P1 ||
    true_donor %in% d_P3 & true_recip %in% a_P2 ||
    true_donor %in% d_P4 & true_recip %in% a_P2 ||
    true_donor %in% d_P4 & true_recip %in% a_P1) {
      results[i,9] <- 2
    } else {
      results[i,9] <- 3
    }
  } else {
    if (test[[8]] < 0.05) {
      if (true_donor %in% d_P1 & true_recip %in% a_P3 ||
      true_donor %in% d_P2 & true_recip %in% a_P3 ||
      true_donor %in% d_P3 & true_recip %in% a_P1 ||
      true_donor %in% d_P3 & true_recip %in% a_P2) {
        results[i,9] <- 4
      } else if (true_donor %in% d_P4 & true_recip %in% a_P2 ||
      true_donor %in% d_P4 & true_recip %in% a_P1){
        results[i,9] <- 5
      } else if (!true_donor %in% d_quatuor & true_recip %in% a_P1 ||
      !true_donor %in% d_quatuor & true_recip %in% a_P2) {
        results[i,9] <- 6
      } else {
        results[i,9] <- 7
      }
    } else {
      results[i,9] <- 1
    }
  }
}



results$H0 <- factor(results$H0, levels = c("0", "1", "2", "3", "4", "5", "6", "7"))
table(results$H0)


all =
relax =
strict =



# is_H_true <- function(test, tree) {
#   a_P1 <- c(test[[1]], Ancestors(tree, test[[1]], type = "all"))
#   a_P2 <- c(test[[2]], Ancestors(tree, test[[2]], type = "all"))
#   a_P3 <- c(test[[3]], Ancestors(tree, test[[3]], type = "all"))
#   first_son_P3 <- a_P3[match(mrca.phylo(tree, c(test[[3]], test[[2]])), a_P3) - 1]
#   if (test[[7]] < 0 & test[[8]] < 0.001) {
#     first_son_P1 <- a_P1[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P1) - 1]
#     if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
#         true_recip %in% a_P1 ||
#         true_donor %in% Descendants(tree, first_son_P1, type = "all") &
#         true_recip %in% a_P3) {return("True")} else {return("False")}
#   } else if (test[[7]] > 0 & test[[8]] < 0.001) {
#     first_son_P2 <- a_P2[match(mrca.phylo(tree, c(test[[1]], test[[2]])), a_P2) - 1]
#     if (true_donor %in% Descendants(tree, first_son_P3, type = "all") &
#         true_recip %in% a_P2 ||
#         true_donor %in% Descendants(tree, first_son_P2, type = "all") &
#         true_recip %in% a_P3) {return("True")} else {return("False")}
#   } else {return("0")}
# }



#
# top = c(named(topologies[1,117]),named(topologies[2,117]),named(topologies[3,117]),named(topologies[4,117]))
# top
#
# spnd[mrca.phylo(tree, c(top[1], top[2]))]
# dist_to_root[mrca.phylo(tree, c(top[1], top[2]))] > dist_to_root[mrca.phylo(tree, c(top[1], top[3]))] & dist_to_root[mrca.phylo(tree, c(top[1], top[3]))] > dist_to_root[mrca.phylo(tree, c(top[1], top[4]))]
# dist_to_root[mrca.phylo(tree, c(top[1], top[3]))] > dist_to_root[mrca.phylo(tree, c(top[1], top[4]))]
# dist_to_root[mrca.phylo(tree, c(top[1], top[4]))]
#
# rem <- c()
#
#
# p1 <- 16
# p2 <- 14
# p3 <- 5
# p4 <- 1
# (d12 <- dist_to_root[mrca.phylo(tree, c(p1, p2))])
# (d13 <- dist_to_root[mrca.phylo(tree, c(p1, p3))])
# (d14 <- dist_to_root[mrca.phylo(tree, c(p1, p4))])
# d12 > d13 & d12 > d14 & d13 > d14
#
# l1=l2=l3=l4 = seq(from = 1, to = nrow(sites))
# comb <- t(CJ(l1,l2,l3,l4)[l1 != l2 & l1 != l3 & l1 != l4 & l2 != l3 & l2 != l4 & l3 != l4])
# rem <- c()
# for (i in 1:ncol(comb)) {
#   p1 <- named(comb[1,i])
#   p2 <- named(comb[2,i])
#   p3 <- named(comb[3,i])
#   p4 <- named(comb[4,i])
#   d12 <- dist_to_root[mrca.phylo(tree, c(p1, p2))]
#   d13 <- dist_to_root[mrca.phylo(tree, c(p1, p3))]
#   d14 <- dist_to_root[mrca.phylo(tree, c(p1, p4))]
#   if (d12 > d13 & d12 > d14 & d13 > d14) {
#     rem <- c(rem, i)
# }}
# comb_f = t(unique(apply(comb[,rem], 2, function(x) {return(c(min(x[c(1,2)]), max(x[c(1,2)]), x[3], x[4]))})))
# (topologies <- t(unique(comb_f[duplicated(comb_f),])))
#
#
# l1=l2=l3=l4 = seq(from = 1, to = nrow(sites))
# comb <- t(CJ(l1,l2,l3,l4)[l1 != l2 & l1 != l3 & l1 != l4 & l2 != l3 & l2 != l4 & l3 != l4])
# rem <- c()
# for (i in 1:ncol(comb)) {
#   if (dist_to_root[mrca.phylo(tree, c(named(comb[1,i]), named(comb[2,i])))] > dist_to_root[mrca.phylo(tree, c(named(comb[2,i]), named(comb[3,i])))] &
#   dist_to_root[mrca.phylo(tree, c(named(comb[2,i]), named(comb[3,i])))] > dist_to_root[mrca.phylo(tree, c(named(comb[3,i]), named(comb[4,i])))]) {
#     rem <- c(rem, i)
# }}
# comb_f = t(unique(apply(comb[,rem], 2, function(x) {return(c(min(x[c(1,2)]), max(x[c(1,2)]), x[3], x[4]))})))
# topologies <- t(unique(comb_f[duplicated(comb_f),]))
