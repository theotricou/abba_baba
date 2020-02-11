#!/usr/bin/env Rscript
# Théo
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (ms command file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "output"
}

require('coala')
require('phyclust')

activate_ms(priority = 600)

is_abba <- function(seg_site) {if (seg_site[1] == seg_site[4] & seg_site[2] == seg_site[3] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

is_baba <- function(seg_site) {if (seg_site[1] == seg_site[3] & seg_site[2] == seg_site[4] & seg_site[1] != seg_site[2]) {return(1)} else {return(0)}}

named <- function(P) {return(which(pop == P)[1])}

D_stat <- function(stat_simulation, quatuor){
  abba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_abba(x)))
  baba <- sum(apply(sites[c(quatuor[1],quatuor[2],quatuor[3],quatuor[4]),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  data <- list("quatuor" = c(named(quatuor[1]), named(quatuor[2]), named(quatuor[3]), named(quatuor[4])), "D-stat" = D, "abba" = abba, "baba" = baba)
  return(data)}

file = args[1]
command <- readChar('ms_command.txt', file.info('ms_command.txt')$size)
command_exe = parse(text = command)
model <- eval(command_exe)
pop <- eval(parse(text = strsplit(strsplit(command, ", loci")[[1]][1], " = ")[[1]][2]))
for (i in 2:length(pop)) {pop[i] <- pop[i-1] + pop[i]}

rep = simulate(model, nsim = 10000, core = 1)

uniq <- lapply(rep, function(x){
  if (ncol(x$seg_sites[[1]][[1]]) == 1) {
      return(x$seg_sites[[1]][[1]])}})

sites <- c()
for (i in 1:length(uniq)) {
  if (length(uniq[[i]]) != 0) {
    if (length(sites) == 0) {
      sites <- cbind(uniq[[i]])
    } else {
      sites <- cbind(sites, uniq[[i]])}}}

comb = combn(seq(from = 1, to = nrow(sites)),4) # retunr all possible quatuor (only leaves)

results = as.data.frame(do.call(rbind, apply(comb, 2, function(x) D_stat(sites, x))))

# now test if basic H1 is true: P3 introgressed in P2 (D > 0 or P1 (D < 1)
