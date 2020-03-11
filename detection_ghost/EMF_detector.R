#!/usr/bin/env Rscript)
# Théo

for i in tes*/data*; do
  sed "1d" $i >> second
done

R


data = read.table('second')
colnames(data) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")

sub = data[data$Pvalue < 0.05,]
options(scipen=999)

min = seq(0,990000 ,100000)
max = seq(100000, 1000000, 100000)
level = as.factor(paste(min, max, sep='\n'))
par(mfrow=c(1,3))

H1 = c()
H2 = c()
H2e = c()
H2c = c()
for (i in 1:length(min)) {
  a = sub[sub$d_p13_p14 > min[i] & sub$d_p13_p14 <= max[i],]
  H1[i] = sum(nrow(a[a$type_D == "P3" & a$type_R == "P1",]),
  nrow(a[a$type_D == "P3" & a$type_R == "P2",]),
  nrow(a[a$type_D == "P1" & a$type_R == "P3",]),
  nrow(a[a$type_D == "P2" & a$type_R == "P3",]))
  H2[i] = sum(nrow(a[a$type_D == "N2" & a$type_R == "P1",]),
  nrow(a[a$type_D == "N2" & a$type_R == "P2",]))
  H2e[i] = sum(H2[i],
  nrow(a[a$type_D == "P4" & a$type_R == "P1",]),
  nrow(a[a$type_D == "P4" & a$type_R == "P2",]))
  H2c[i] = sum(H2e[i],
  nrow(a[a$type_D == "N3" & a$type_R == "P1",]),
  nrow(a[a$type_D == "N3" & a$type_R == "P2",]))
}
per_error=H2/(H1+H2)
per_error_ext=H2e/(H1+H2e)
per_error_com=H2c/(H1+H2c)
rate <- rbind(per_error, per_error_ext, per_error_com)
mp <- barplot(rate, beside = TRUE, cex = 0.75, cex.lab  = 1 , cex.axis = 0.7, las = 1,
              ylim = c(0, 1.1),
              main="ABBA/BABA test error rate",
              names.arg = level,
              xlab = "Distance between P3 and P4",
              ylab = "Error rata",
              col = c("red", "blue", "green"))
legend("topleft",
       legend = c("H2 error rate","H2ext error rate","H2com error rate"),
       fill = c("red", "blue", "green"))





### diversification comparaison

data00 = read.table('00_div/second')
data01 = read.table('01_div/second')
data05 = read.table('05_div/second')
data09 = read.table('09_div/second')
colnames(data00) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")
colnames(data01) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")
colnames(data05) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")
colnames(data09) <- c("P1","P2","P3","P4","abba","baba","Dstat","Pvalue","d_p13_p14","type_D","type_R")
options(scipen=999)

sub00 = data00[data00$Pvalue < 0.05,]
sub01 = data01[data01$Pvalue < 0.05,]
sub05 = data05[data05$Pvalue < 0.05,]
sub09 = data09[data09$Pvalue < 0.05,]

bin = 200000

min = seq(0,1000000 - bin ,bin)
max = seq(bin, 1000000, bin)
level = as.factor(paste(min, max, sep='\n'))

H100 =  c()
H200 =  c()
H101 =  c()
H201 =  c()
H205 =  c()
H105 =  c()
H109 = c()
H209 = c()
for (i in 1:length(min)) {
  a00 = sub00[sub00$d_p13_p14 > min[i] & sub00$d_p13_p14 <= max[i],]
  H100[i] = sum(nrow(a00[a00$type_D == "P3" & a00$type_R == "P1",]),
  nrow(a00[a00$type_D == "P3" & a00$type_R == "P2",]),
  nrow(a00[a00$type_D == "P1" & a00$type_R == "P3",]),
  nrow(a00[a00$type_D == "P2" & a00$type_R == "P3",]))
  H200[i] = sum(nrow(a00[a00$type_D == "N2" & a00$type_R == "P1",]),
  nrow(a00[a00$type_D == "N2" & a00$type_R == "P2",]))

  a01 = sub01[sub01$d_p13_p14 > min[i] & sub01$d_p13_p14 <= max[i],]
  H101[i] = sum(nrow(a01[a01$type_D == "P3" & a01$type_R == "P1",]),
  nrow(a01[a01$type_D == "P3" & a01$type_R == "P2",]),
  nrow(a01[a01$type_D == "P1" & a01$type_R == "P3",]),
  nrow(a01[a01$type_D == "P2" & a01$type_R == "P3",]))
  H201[i] = sum(nrow(a01[a01$type_D == "N2" & a01$type_R == "P1",]),
  nrow(a01[a01$type_D == "N2" & a01$type_R == "P2",]))

  a05 = sub05[sub05$d_p13_p14 > min[i] & sub05$d_p13_p14 <= max[i],]
  H105[i] = sum(nrow(a05[a05$type_D == "P3" & a05$type_R == "P1",]),
  nrow(a05[a05$type_D == "P3" & a05$type_R == "P2",]),
  nrow(a05[a05$type_D == "P1" & a05$type_R == "P3",]),
  nrow(a05[a05$type_D == "P2" & a05$type_R == "P3",]))
  H205[i] = sum(nrow(a05[a05$type_D == "N2" & a05$type_R == "P1",]),
  nrow(a05[a05$type_D == "N2" & a05$type_R == "P2",]))

  a09 = sub09[sub09$d_p13_p14 > min[i] & sub09$d_p13_p14 <= max[i],]
  H109[i] = sum(nrow(a09[a09$type_D == "P3" & a09$type_R == "P1",]),
  nrow(a09[a09$type_D == "P3" & a09$type_R == "P2",]),
  nrow(a09[a09$type_D == "P1" & a09$type_R == "P3",]),
  nrow(a09[a09$type_D == "P2" & a09$type_R == "P3",]))
  H209[i] = sum(nrow(a09[a09$type_D == "N2" & a09$type_R == "P1",]),
  nrow(a09[a09$type_D == "N2" & a09$type_R == "P2",]))

}

per_error00=H200/(H100+H200)
per_error01=H201/(H101+H201)
per_error05=H205/(H105+H205)
per_error09=H209/(H109+H209)


rate <- rbind(per_error00, per_error01, per_error05, per_error09)

mp <- barplot(rate, beside = TRUE, cex = 0.75, cex.lab  = 1 , cex.axis = 0.7, las = 1,
              ylim = c(0, 1.1),
              main="ABBA/BABA test error rate",
              names.arg = level,
              xlab = "Distance between P3 and P4",
              ylab = "Error rata",
              col = c("red", "blue", "green", "black"))
legend("topleft", title = "Diversification rate (extinction rate/speciation rate)",
       legend = c("0,0","0,1","0,5", "0,9"),
       fill = c("red", "blue", "green", "black"))








df <- data.frame(level=factor(level),outside=per_error)
plot(df, col = "red", cex = 0.1, cex.lab  = 1 , cex.axis = 0.7, las = 1,
xlab = "Distance between P3 and P4",
ylab = "Error rata",
main = "ABBA/BABA test error rate")
dfe <- data.frame(level=factor(level),outside=per_error_ext)
plot(dfe, col = "red", cex = 0.1, cex.lab  = 1 , cex.axis = 0.7, las = 1,
xlab = "Distance between P3 and P4",
ylab = "Error rata",
main = "ABBA/BABA test error rate")
dfc <- data.frame(level=factor(level),outside=per_error_com)
plot(df, col = "red", cex = 0.1, cex.lab  = 1 , cex.axis = 0.7, las = 1,
xlab = "Distance between P3 and P4",
ylab = "Error rata",
main = "ABBA/BABA test error rate")



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




sum(nrow(sub[sub$type_D == "N2" & sub$type_R == "P1",]),
  nrow(sub[sub$type_D == "N2" & sub$type_R == "P2",]))



sum(nrow(sub[sub$type_D == "P3" & sub$type_R == "P1",]),
  nrow(sub[sub$type_D == "P3" & sub$type_R == "P2",]),
  nrow(sub[sub$type_D == "P1" & sub$type_R == "P3",]),
  nrow(sub[sub$type_D == "P2" & sub$type_R == "P3",]))



library('ggplot2')


data = as.data.frame(t(rbind(levels(level), per_error)))
colnames(data) = c("range","error")
data$error = as.numeric(as.character(data$error))

ggplot(data=data, aes(x=range, y=error)) + geom_boxplot() +
coord_cartesian(ylim = c(0, 1))

ggplot(data=data,                             # Define the data to plot
       aes(x=vore,
           fill=vore)) +                        # Define how bars are coloured
  geom_bar() +
  scale_x_discrete(name = 'Feeding Type',
                   labels=labs) +
  scale_fill_discrete(name = 'Feeding Type',
                      labels=labs) +
  labs(x='Feeding Type',
       y='Number of Species') +
  theme_bw()

Ordonnons les modalités d'absentéisme
etudiants$abs <- ordered(etudiants$abs, levels = c("0-1", ">1"))
par(mfrow = c(1, 2))
# Pour avoir deux graphiques côte à côte
boxplot(etudiants$note~etudiants$mBac+etudiants$abs, main = "Mention au bac et absentéisme",
xlab = "Mention au bac et absentéisme", ylab = "Note en MathSV sur 20", las = 1,
 col = c("pink", "grey", "lightblue"), varwidth = TRUE, notch = TRUE)
 boxplot(etudiants$note~etudiants$abs+etudiants$mBac, main = "Absentéisme et mention au bac",
  xlab = "Absentéisme et mention au bac", ylab = "Note en MathSV sur 20", las = 1,
   col = rep(c("pink", "grey", "lightblue"), each = 2), varwidth = TRUE, notch = TRUE)


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
