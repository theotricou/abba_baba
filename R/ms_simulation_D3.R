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

devil_mambojambo<-function(range, trees){
  return(Reduce('+', lapply(trees[as.integer(fromto[range,1]):as.integer(fromto[range,2])], function(x) cophenetic(x)[ord,ord])))
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
  if (true_donor %in% leaf_descendants(p3,p2)){
    return("P3")
  }else if (true_donor %in% leaf_descendants(p2,p1)){
    return("P2")
  }else if (true_donor %in% leaf_descendants(p1,p2)){
    return("P1")
  }else if (true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){
    return("N1")
  }else{
    return("O")
  }
}

position_R_in_tree<-function(p1,p2,p3){
  if (true_recip %in% leaf_ancestors(p3,p2)){
    return("P3")
  }else if (true_recip %in% leaf_ancestors(p2,p1)){
    return("P2")
  }else if (true_recip %in% leaf_ancestors(p1,p2)){
    return("P1")
  }else if (true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){
    return("N1")
  }else{
    return("O")
  }
}

D3_v2<-function(matrice, t){
  P1<-da[da$nu == t[1],3]
  P2<-da[da$nu == t[2],3]
  P3<-da[da$nu == t[3],3]
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  D3_donor<-position_D_in_tree(t[1],t[2],t[3])
  D3_recip<-position_R_in_tree(t[1],t[2],t[3])
  # sister?
  donor<-D3_donor
  temp<-c(t[1],t[2],t[3])
  if (donor %in% c("O",'N2','N1')){
    sista = 'Nope'
  }else{
    col<-as.numeric(strsplit(donor,'P')[[1]][2])
    node=as.numeric(temp[col])
    if (true_donor %in% c(node, Ancestors(tree, node))){
      sista = 'True'
    }else{
      sista = 'Sister'
    }
  }
  desc<-spnd[Descendants(tree, mrca.phylo(tree, c(t[1], t[3])))[[1]]]
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
    "dP1P3" = as.integer(strsplit(spnd[mrca.phylo(tree, c(t[1], t[3]))], "@")[[1]][2]),
    "SisDonor" = sista,
    "ali"=length(grep("ali>", desc)),
    "uns"=length(grep("uns>", desc)),
    "ext"=length(grep("ext>", desc))
  )
  return(data)
}

D3_boot<-function(matrice, t){
  P1<-da[da$nu == t[1],3]
  P2<-da[da$nu == t[2],3]
  P3<-da[da$nu == t[3],3]
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  return(as.numeric(D3))
}

bootstraps<-function(trees){
    boot_trees<-sample(trees, 1000, replace=TRUE)
    blocks<-seq(1,length(boot_trees), by=50000)
    fromto<-cbind(blocks,c(blocks[-1]-1, length(boot_trees)))
    devil_mambojambo<-function(range, trees){
      return(Reduce('+', lapply(trees[as.integer(fromto[range,1]):as.integer(fromto[range,2])], function(x) cophenetic(x)[ord,ord])))
    }
    list_submat<-lapply(1:length(blocks),function(x) devil_mambojambo(x, boot_trees))
    sum_mat<-Reduce('+',list_submat)

    return(do.call(rbind,mclapply(as.data.frame(t(topologiesD3)), function(x) D3_boot(sum_mat, x),mc.cores = N_CORE)))
}

Zscore<-function(D3,boot_D3){
  D_Z      <- D3 / sd(boot_D3, na.rm=T)
  return(D_Z)
}

################################################################################

### D3 stat

tree <- read.tree(file.path(output, "spe_tree"))
extant<-grep("ali>",tree$tip.label)
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

spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

single_trees<-lapply(unlist(lapply(which(rep != "NULL"), function(x) rep[[x]]$trees[[1]])), function(x) read.tree(text=x))

topologiesD3 <- gettrios(tr)
topologiesD3 <- apply(topologiesD3,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
colnames(topologiesD3)<-c('P1','P2','P3')

ord<-sort(row.names(cophenetic(single_trees[[1]])))
cat("step 3 : computing D3\n")

blocks<-seq(1,length(single_trees), by=50000)
fromto<-cbind(blocks,c(blocks[-1]-1, length(single_trees)))
list_submat<-lapply(1:length(blocks),function(x) devil_mambojambo(x, single_trees))
sum_mat<-Reduce('+',list_submat)

results<-as.data.frame(do.call(rbind,mclapply(as.data.frame(t(topologiesD3)), function(x) D3_v2(sum_mat, x),mc.cores = N_CORE)))
results$allsp=length(tree$tip.label)

cat("step 4 : Bootstraps D3\n")

boot=do.call('cbind',mclapply(1:1000, function(x) bootstraps(single_trees),mc.cores = N_CORE))
results$Z_score<-unlist(lapply(1:nrow(results), function(x) Zscore(as.numeric(as.character(results[x,"D3"])), unlist(boot[x,]))))

cat("step 5 : output D3\n")

outfile <- paste(output, "data_D3.txt", sep = "/")
write.table(results, outfile, sep = "\t", row.names = F, append = F, quote=F)

cat("End \n")
