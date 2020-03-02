#!/usr/bin/env Rscript)
# Théo


data = read.table('data.txt', h = T)
d = data[data$Pvalue < 0.001,]

P_lins <- unique(c(
paste("N",levels(as.factor(data$P1)), sep=""),
paste("N",levels(as.factor(data$P2)), sep=""),
paste("N",levels(as.factor(data$P3)), sep=""),
paste("N",levels(as.factor(data$P4)), sep="")
))

N_lins <- paste("N",unique(paste(data$P3, data$P4, spe = "_"), sep="_")
O_lins <- paste("N",unique(data$P4), "O", sep = "_")
lineages <- sort(c(P_lins, N_lins, O_lins))


matrix <- matrix(0, nrow = length(lineages), ncol = length(lineages))
colnames(matrix) <- lineages
rownames(matrix) <- lineages


for (i in 1:nrow(d)) {
  line = d[i,]
  pp = paste("N", line[1:4], sep = "")
  nn = c(paste("N", pp$P3, pp$P4 , sep = "_"), paste("N", line$P4, "O", sep = "_") )

  paste("N", c((line$P3 * 10000) + line$P4 , (line$P4 * 100)), sep = "")
  if (line[7] < 0) {
    pp = c(pp[2],pp[1],pp[3],pp[4])
  }
  matrix[pp[4],pp[1]] = matrix[pp[4],pp[1]] + 1
  matrix[pp[1],pp[4]] = matrix[pp[1],pp[4]] + 1
  matrix[pp[3],pp[2]] = matrix[pp[3],pp[2]] + 1
  matrix[pp[2],pp[3]] = matrix[pp[2],pp[3]] + 1
  matrix[nn[1],pp[1]] = matrix[nn[1],pp[1]] + 1
  matrix[nn[2],pp[1]] = matrix[nn[2],pp[1]] + 1
}


sum(matrix)
max(matrix)

temp = matrix
temp = temp[rowSums(temp[,-1]) != 0,]
temp = temp[,colSums(temp[-1,]) != 0]

# heatmap(t(temp), Colv = NA, Rowv = NA, scale="column", col= colorRampPalette(brewer.pal(8, "Blues"))(25))

library(RColorBrewer)
col<- colorRampPalette(c("white", "lightgrey ", "red"))(256)
heatmap(data, Colv = NA, Rowv = NA, scale="none", col = col)
