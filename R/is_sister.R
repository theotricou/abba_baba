#!/usr/bin/env Rscript
# Theo
# donor is sister ?

library(ape)
library(phangorn)

source2 <- function(file, start, end, ...) {
    file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
    file.lines.collapsed <- paste(file.lines, collapse='\n')
    source(textConnection(file.lines.collapsed), ...)
}

source2('ms_command.R', 7,10 )

results<-read.table('data.txt', h=T)

tree <- read.tree("spe_tree")
extant<-grep("l>",tree$tip.label)
tr<-keep.tip(tree, tree$tip.label[extant])
da<-data.frame(na=tree$tip.label[extant], nu=extant, ge=paste("s",1:length(extant),sep=""))
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file


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

results$is_sister<-apply(results,1, function(x) is_sister(x))
write.table(results, "data_s.txt", sep = "\t", row.names = F, append = F, quote=F)
