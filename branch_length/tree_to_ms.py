#! /usr/bin/env python3
# Theo Tricou

import os
import sys
from ete3 import Tree as tr
import random
import argparse
import numpy as np

def gene_t(node):
    return(str(node.name).split('@')[1])

def gene_n(node):
    if ">" in node.name:
        return(str(node.name).split('@')[0].split('>')[1])
    else:
        return(str(node.name).split('@')[0])

def tree_new_dist(tree, n_generation_to_root):
    multi = n_generation_to_root / tree.get_farthest_leaf()[1]
    for i in t.traverse():
        i.dist = int((i.dist * multi)//1)
    return(tree)

t = tr("tree", format = 1) # read phylo tree
sampled = ["90093", "90106", "90110"]

Ne = 100_000
len_locus = 1
ploidy = 1
n_generation_to_root = 1_000_000

t = tree_new_dist(t, n_generation_to_root)

Head = D_R = Pop = Coal = Merge = Migration_starts = Migration_ends = Stat_sum = True_migration = str() # coal_model extended parameters
ext_lineages = []
count = 1
trio= []
for i in t.traverse('postorder'):
    old=i.name
    if i.is_leaf():
        if i.name in sampled:
            ext_lineages.append(1)
            i.name = "ali>" + str(count) + "@0"
        else:
            ext_lineages.append(0)
            i.name = "ext>" + str(count) + "@" + str(int(n_generation_to_root - i.get_distance(t)))
        count += 1
    else:
        pop_d = str(gene_n(i.get_descendants()[0])).split("_")[0]
        pop_r = str(gene_n(i.get_descendants()[1])).split("_")[0]
        pop_g = int(n_generation_to_root - i.get_distance(t))
        i.name = pop_d + "_" + pop_r + "@" + str(pop_g)
        pop_t =  pop_g / (4 * Ne)
        Merge += "feat_pop_merge(%s, %s, %s) + " % (str(pop_t), pop_r, pop_d) + "\n"
    if old in sampled:
        trio.append(i.name)

t.write(outfile = 'spe_tree', format=1, format_root_node=True)

# P3 , P2, P1
# 'ali>25@0', 'ali>32@0', 'ali>36@0'

migration_fraction = 0.2
migration_rate = migration_fraction / migration_time
migration_start = 0
migration_end = migration_start + migration_time

# R source code variables
Head = "library('ape') \n" + "library('coala') \n" + "library('phyclust') \n" + "library('phangorn') \n" + "activate_ms(priority = 600) \n\n"
Coal = "model <- coal_model(sample_size = c(%s), loci_number = %s, loci_length = 1, ploidy = %s) + \n" % (', '.join(str(x) for x in ext_lineages), len_locus, ploidy)
Mutation = "feat_mutation(rate = 20.0, model = 'IFS', fixed_number = FALSE, locus_group = 'all') + \n"
Stat_sum = "sumstat_seg_sites() + sumstat_trees() \n\n"


# Ingroup introgression between P2 and P3
pop_donor = t&"ali>32@0"
pop_recip = t&'ali>25@0'
ingroup_introgression = """
model_ingroup<-model +
feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') +
feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') \n
""" % (migration_rate, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0],
migration_start, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0],
migration_end)


# Ghost introgression between a ghost and P2
pop_recip = t&"ali>32@0"
pop_donor = t&'ext>52@11'
ghost_introgression = """
model_ghost<-model +
feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') +
feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') \n
""" % (migration_rate, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0],
migration_start, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0],
migration_end)


with open('ms_command.R', "w") as output:
    output.write("".join([Head, Coal, Mutation, Merge, Stat_sum, ingroup_introgression, ghost_introgression]))

#GNU Terry Pratchett
