#!/usr/bin/env Rscript
# Théo


cat("Starting ms simulation")
args = commandArgs(trailingOnly=TRUE)

cat("step 1 : initializing \n")

# source2 <- function(file, start, end, ...) {
#     file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
#     file.lines.collapsed <- paste(file.lines, collapse='\n')
#     source(textConnection(file.lines.collapsed), ...)
# }

# if (length(args) == 0) {
#   if(file.exists("Simulation/ms_command.R")) {
#     source("Simulation/ms_command.R")
#     source("Parameters/sim_parameters")
#     output = "Simulation"
#   } else {
#     print("You need to create a ms command by running build_ms_command.py to run this script")
#     quit()
#   }
# } else if (length(args) == 1) {
#   source(args[1])
#   path = strsplit(args[1], "/")[[1]]
#   output = paste(path[seq(1,length(path) -1)], collapse = '/')
#   source2(paste(output, paste("Param_", output, sep = ""), sep = "/"),37,85)
#   if (length(path) == 1) {output = "./"} else {output = paste(path[seq(1,length(path) -1)], collapse = '/')}
# } else if (length(args) == 2) {
#   source(args[1])
#   output = args[2]
# } else {
#   print("Too much argument were given, only the two first one will be used as input and output respectively")
#   source(args[1])
#   path = strsplit(args[2], "/")[[1]]
#   if (length(path) == 1) {output = "./"} else {output = paste(path[seq(1,length(path) -1)], collapse = '/')}
# }

# source(args[1])
# source(args[2])
# path = strsplit(args[1], "/")[[1]]
# output = paste(path[seq(1,length(path) -1)], collapse = '/')

output=args[1]
source(paste(args[1],"ms_command.R", sep ="/"))
source(paste(args[1],"sim_parameters", sep ="/"))


is_abba <- function(seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[4]){return(1)} else {return(0)}}
is_baba <- function(seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[3]){return(1)} else {return(0)}}


validateandreorder<-function(arr, dist) {
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  } else if (sum(submat==max(submat))!=6 & sum(submat==max(submat))!=8) {
    return("ERROR in rounding branch length")
  }
}

getquatuors<-function(tr) {
  dist<-cophenetic(tr)
  dist<-round(dist,4)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}

leafs_relatives <- function(tree, leaf, comp_leaf) {
  all_anc_leaf <- c(leaf, Ancestors(tree, leaf, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(tree, c(leaf, comp_leaf)), all_anc_leaf) - 1]
  # only_anc_leaf <- all_anc_leaf[1:match(first_son_leaf, all_anc_leaf)]
  des_leaf <- Descendants(tree, first_son_leaf, type = "all")
  result <- list("a_anc" = all_anc_leaf, "firs" = first_son_leaf,  "des" = des_leaf) #"o_anc" = only_anc_leaf,
  return(result)
}

where_is_Daldo <- function(tree, donor, p1, p2, p3, p4) {
  if (donor == p1) { return("P1")}
  else if (donor == p2) { return("P2")}
  else if (donor == p3) { return("P3")}
  else if (donor == p4) { return("P4")}
  else {
    if (!donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")) {return("N3")}
    else {
      if (!donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")) {
        if (donor %in% leafs_relatives(tree, p4, p1)$des) {return("P4")}
        else {return("N2")}
      } else {
        if (!donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p2)), type="all")) {
          if (donor %in% leafs_relatives(tree, p3, p1)$des) {return("P3")}
          else {return("N1")}
        } else {
          if (donor %in% leafs_relatives(tree, p2, p1)$des) {return("P2")}
          else {return("P1")}
}}}}}

where_is_Raldo <- function(tree, recip, p1, p2, p3, p4) {
  if (recip == p1) { return("P1")}
  else if (recip == p2) { return("P2")}
  else if (recip == p3) { return("P3")}
  else if (recip == p4) { return("P4")}
  else {
    if (!recip %in% Descendants(tree, mrca.phylo(tree, c(p1, p4)), type="all")){ return("N3")}
    else {
      if (!recip %in% Ancestors(tree, p4) &
      !recip %in% Ancestors(tree, p3) &
      !recip %in% Ancestors(tree, p2) &
      !recip %in% Ancestors(tree, p1)) {return("N3")}
      else{
        if (!recip %in% Descendants(tree, mrca.phylo(tree, c(p1, p3)), type="all")) {
          if (recip %in% c(p4, Ancestors(tree, p4))) {return("P4")}
          else {return("N2")}
        } else {
          if (!recip %in% Descendants(tree, mrca.phylo(tree, c(p1, p2)), type="all")) {
            if (recip %in% c(p3, Ancestors(tree, p3))) {return("P3")}
            else {return("N1")}
          } else {
            if (recip %in% c(p1, Ancestors(tree, p1))) {return("P1")}
            else if (recip %in% c(p2, Ancestors(tree, p2))) {return("P2")}
}}}}}}

D_stat <- function(stat_simulation, quatuor){
  P1 = as.integer(strsplit(quatuor[1], "@")[[1]][1])
  p1 = pop[P1]
  P2 = as.integer(strsplit(quatuor[2], "@")[[1]][1])
  p2 = pop[P2]
  P3 = as.integer(strsplit(quatuor[3], "@")[[1]][1])
  p3 = pop[P3]
  P4 = as.integer(strsplit(quatuor[4], "@")[[1]][1])
  p4 = pop[P4]
  abba <- sum(apply(sites[c(p1, p2, p3, p4),], 2, function(x) is_abba(x)))
  baba <- sum(apply(sites[c(p1, p2, p3, p4),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  ng_p3 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(P1, P3))], "@")[[1]][2])
  ng_p4 <- as.integer(strsplit(spnd[mrca.phylo(tree, c(P1, P4))], "@")[[1]][2])
  dist_p13_p14 <- ng_p4 - ng_p3
  donor = where_is_Daldo(tree, true_donor, P1, P2, P3, P4)
  recip = where_is_Raldo(tree, true_recip, P1, P2, P3, P4)
  data <- c("P1" = P1, "P2" = P2, "P3" = P3, "P4" = P4,
  "abba" = abba, "baba" = baba, "Dstat" = D,
  "Pvalue" = binom.test(c(abba, baba), p = 0.5, conf.level= 0.99)$p.value,
  "d_p13_p14" = dist_p13_p14,
  "type_D" = donor, "type_R" = recip)
  return(data)
}

# output = "sim_1"
# source(paste(output, "ms_command.R", sep = "/"))
# source(paste(output, "Parameters", sep = "/"))

tree <- read.tree(file.path(output, "spe_tree"))
tr = keep.tip(tree, tree$tip.label[which(pop == 1)])

for (i in 2:length(pop)) {pop[i] <- pop[i-1] + pop[i]}

cat("step 2 : running simmulation \n")
# rep = simulate(model, nsim = 100000, cores = 4, seed = 23805)

if (SEED == 0) {
  rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
} else {
  rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE, seed = SEED)
}



coal_trees <- c()
for (i in 1:length(rep)) {coal_trees[i] <- rep[[i]]$trees[[1]]}
cat("Outputting all trees from simualtion in: ")
outfile_a <- paste(output, "all_trees", sep = "/")
write(coal_trees, file=outfile_a)
cmd <- paste("tar -cvzf", paste(outfile_a,"tar.gz", sep = "."), outfile_a, "--remove-files", sep = " ")
system(cmd, wait=T)

uniq <- lapply(rep, function(x){
  if (ncol(x$seg_sites[[1]][[1]]) == 1) {
      return(x$seg_sites[[1]][[1]])}})

cat("Outputting trees from single segregating site locus in: ")
outfile_s <- paste(output, "single_trees", sep = "/")
single_trees <- coal_trees[which(uniq != "NULL")]
write(single_trees, file=outfile_s)
cmd <- paste("tar -cvzf", paste(outfile_s,"tar.gz", sep = "."), outfile_s, "--remove-files", sep = " ")
system(cmd, wait=T)


sites <- c()
for (i in 1:length(uniq)) {
  if (length(uniq[[i]]) != 0 & sum(uniq[[i]]) > 1) {
    sites <- cbind(sites, uniq[[i]])}}

topologies <- getquatuors(tr)
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

cat("step 3 : computing summary statistics \n")

results <- as.data.frame(t(apply(topologies, 1, function(x) D_stat(sites, x)))) # bottleneck

cat("step 4 : output \n")

outfile_d <- paste(output, "data.txt", sep = "/")
write.table(results, outfile_d, sep = "\t", row.names = F, append = F, quote=F)

cat("End \n")
