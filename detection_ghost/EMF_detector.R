#!/usr/bin/env Rscript)
# Théo

# global stat
data = read.table('second')
colnames(data) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")

a = data[data$Pvalue < 0.05,]
# true positives H1
nrow(a[a$type_D == "P3" & a$type_R == "P1",])
nrow(a[a$type_D == "P3" & a$type_R == "P2",])
nrow(a[a$type_D == "P1" & a$type_R == "P3",])
nrow(a[a$type_D == "P2" & a$type_R == "P3",])

sum(nrow(a[a$type_D == "P3" & a$type_R == "P1",]),
nrow(a[a$type_D == "P3" & a$type_R == "P2",]),
nrow(a[a$type_D == "P1" & a$type_R == "P3",]),
nrow(a[a$type_D == "P2" & a$type_R == "P3",]))

# true positives H2 strict
nrow(a[a$type_D == "N2" & a$type_R == "P1",])
nrow(a[a$type_D == "N2" & a$type_R == "P2",])

sum(nrow(a[a$type_D == "N2" & a$type_R == "P1",]),
nrow(a[a$type_D == "N2" & a$type_R == "P2",]))

# true positives H2 extanded
nrow(a[a$type_D == "P4" & a$type_R == "P1",])
nrow(a[a$type_D == "P4" & a$type_R == "P2",])

sum(nrow(a[a$type_D == "N2" & a$type_R == "P1",]),
nrow(a[a$type_D == "N2" & a$type_R == "P2",]),
nrow(a[a$type_D == "P4" & a$type_R == "P1",]),
nrow(a[a$type_D == "P4" & a$type_R == "P2",]))

# true positives H2 complete
nrow(a[a$type_D == "N3" & a$type_R == "P1",])
nrow(a[a$type_D == "N3" & a$type_R == "P2",])

sum(nrow(a[a$type_D == "N2" & a$type_R == "P1",]),
nrow(a[a$type_D == "N2" & a$type_R == "P2",]),
nrow(a[a$type_D == "P4" & a$type_R == "P1",]),
nrow(a[a$type_D == "P4" & a$type_R == "P2",]),
nrow(a[a$type_D == "N3" & a$type_R == "P1",]),
nrow(a[a$type_D == "N3" & a$type_R == "P2",]))

nrow(a[a$type_R == "N1",])
nrow(a[a$type_R == "N2",])
nrow(a[a$type_R == "N3",])
nrow(a[a$type_D == "N1",])



unit = sort(unique(data$d_p13_p14))

df <- data.frame(unit=integer(),
vph1=integer(),vph2=integer(),vph2e=integer(),vph2c=integer(),
fph1=integer(),fph2=integer(),fph2e=integer(),fph2c=integer(),
fnh1=integer(),fnh2=integer(),fnh2e=integer(),fnh2c=integer(),
vnh1=integer(),vnh2=integer(),vnh2e=integer(),vnh2c=integer()
)

for (i in 1:length(unit)) {
  tempp = data[data$d_p13_p14 <= unit[i] & data$Pvalue <= 0.001,]
  tempn = data[data$d_p13_p14 <= unit[i] & data$Pvalue > 0.001,]

  df[i, "unit"] = unit[i]

  df[i, "vph1"] = sum(nrow(tempp[tempp$type_D == "P3" & tempp$type_R == "P1",]),
  nrow(tempp[tempp$type_D == "P3" & tempp$type_R == "P2",]),
  nrow(tempp[tempp$type_D == "P1" & tempp$type_R == "P3",]),
  nrow(tempp[tempp$type_D == "P2" & tempp$type_R == "P3",]))
  df[i, "vph2"] = sum(nrow(tempp[tempp$type_D == "N2" & tempp$type_R == "P1",]),
  nrow(tempp[tempp$type_D == "N2" & tempp$type_R == "P2",]))
  df[i, "vph2e"] = sum(df[i, "vph2"], nrow(tempp[tempp$type_D == "P4" & tempp$type_R == "P1",]),
  nrow(tempp[tempp$type_D == "P4" & tempp$type_R == "P2",]))
  df[i, "vph2c"] = sum(df[i, "vph2e"], nrow(tempp[tempp$type_D == "N3" & tempp$type_R == "P1",]),
  nrow(tempp[tempp$type_D == "N3" & tempp$type_R == "P2",]))
  df[i, "fph1"] = nrow(tempp) - df[i, "vph1"]
  df[i, "fph2"] = df[i, "fph1"] - df[i, "vph2"]
  df[i, "fph2e"] = df[i, "fph1"] - df[i, "vph2e"]
  df[i, "fph2c"] = df[i, "fph1"] - df[i, "vph2c"]

  df[i, "fnh1"] = sum(nrow(tempn[tempn$type_D == "P3" & tempn$type_R == "P1",]),
  nrow(tempn[tempn$type_D == "P3" & tempn$type_R == "P2",]),
  nrow(tempn[tempn$type_D == "P1" & tempn$type_R == "P3",]),
  nrow(tempn[tempn$type_D == "P2" & tempn$type_R == "P3",]))
  df[i, "fnh2"] = sum(nrow(tempn[tempn$type_D == "N2" & tempn$type_R == "P1",]),
  nrow(tempn[tempn$type_D == "N2" & tempn$type_R == "P2",]))
  df[i, "fnh2e"] = sum(df[i, "fnh2"], nrow(tempn[tempn$type_D == "P4" & tempn$type_R == "P1",]),
  nrow(tempn[tempn$type_D == "P4" & tempn$type_R == "P2",]))
  df[i, "fnh2c"] = sum(df[i, "fnh2e"], nrow(tempn[tempn$type_D == "N3" & tempn$type_R == "P1",]),
  nrow(tempn[tempn$type_D == "N3" & tempn$type_R == "P2",]))
  df[i, "vnh1"] = nrow(tempn) - df[i, "fnh1"]
  df[i, "vnh2"] = df[i, "vnh1"] - df[i, "fnh2"]
  df[i, "vnh2e"] = df[i, "vnh1"] - df[i, "fnh2e"]
  df[i, "vnh2c"] = df[i, "vnh1"] - df[i, "fnh2c"]
}


head(df,10)

plot(df$unit, (df$vph1), ylim = c(1, (max(df[,c(2:5)]))), pch = 20, cex = 1)
points(df$unit, (df$vph2), col = "red", pch = 20, cex = 1)
points(df$unit, (df$vph2e), col = "blue", pch = 20, cex = 1)
points(df$unit, (df$vph2c), col = "darkgreen", pch = 20, cex = 1)

plot(df$unit, (df$fph1), ylim = c(1, (max(df[,c(6:9)]))), col = "grey", pch = 20, cex = 0.5)
points(df$unit, (df$fph2), col = "orange", pch = 20, cex = 0.5)
points(df$unit, (df$fph2e), col = "lightblue", pch = 20, cex = 0.5)
points(df$unit, (df$fph2c), col = "lightgreen", pch = 20, cex = 0.5)



# all
# vp
plot(df$unit, (df$vph1), ylim = c(1, (max(df[,c(2:13)]))), pch = 20, cex = 0.5)
points(df$unit, (df$vph2), col = "red", pch = 20, cex = 0.5)
points(df$unit, (df$vph2e), col = "blue", pch = 20, cex = 0.5)
points(df$unit, (df$vph2c), col = "darkgreen", pch = 20, cex = 0.5)

# fp
points(df$unit, (df$fph1), col = "grey", pch = 20, cex = 0.5)
points(df$unit, (df$fph2), col = "orange", pch = 20, cex = 0.5)
points(df$unit, (df$fph2e), col = "lightblue", pch = 20, cex = 0.5)
points(df$unit, (df$fph2c), col = "lightgreen", pch = 20, cex = 0.5)

# fn
points(df$unit, (df$fnh1), ylim = c(1, (max(df[,2:5]))), col = "brown", pch = 20, cex = 0.5)
points(df$unit, (df$fnh2), col = "pink", pch = 20, cex = 0.5)
points(df$unit, (df$fnh2e), col = "purple", pch = 20, cex = 0.5)
points(df$unit, (df$fnh2c), col = "yellow", pch = 20, cex = 0.5)

# vn
# points(df$unit, (df$vnh1), ylim = c(1, (max(df[,2:5]))), col = "brown")
# points(df$unit, (df$vnh2), col = "pink")
# points(df$unit, (df$vnh2e), col = "purple")
# points(df$unit, (df$vnh2c), col = "yellow")



points(unit, posH2s, col = "red")
points(unit, posH2e, col = "blue")
points(unit, posH2c, col = "green")
points(unit, false, col = "grey")










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
