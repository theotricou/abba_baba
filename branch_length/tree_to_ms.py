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
        if i.is_leaf() and i.name in extant:
            i.dist = n_generation_to_root - i.up.get_distance(t)
        else :
            i.dist = int((i.dist * multi)//1)
    return(tree)

def ghsot_alive_at_time(tree, recip):
    results = []
    for i in tree.iter_descendants():
        if i.up.get_distance(t) < recip.get_distance(t) and recip.up.get_distance(t) < i.get_distance(t):
            if not "ali" in i.get_descendants():
                results.append(i)
    return(results)

t = tr("test/T/CompleteTree.nwk", format = 1) # read phylo tree
te = tr("test/T/ExtantTree.nwk", format = 1) # read phylo tree

extant = []
for i in te: extant.append(i.name)

sampled = ["90417", "90447", "90448"]

Ne = 100000
mu = 1e-7
len_locus = 1
ploidy = 1
n_generation_to_root = 10_000_000

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
        elif i.name in extant:
            ext_lineages.append(0)
            i.name = "uns>" + str(count) + "@0"
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

pop_recip = t&"ali>138@0"
pop_doror=False
all_donor = ghsot_alive_at_time(t, pop_recip)
pop_donor=random.sample(all_donor,1)[0]

tea_time = int((max(int(pop_recip.name.split("@")[1]), int(pop_donor.name.split("@")[1]))+
    min(int(pop_recip.up.name.split("@")[1]), int(pop_donor.up.name.split("@")[1])))/2)

# migration Parameters
migration_fraction = 0.1 # fraction of the donor population that migrate in recipient population, during migration_generation generations
migration_time = 1 / (4 * Ne) # time of the migration in 4Ne
migration_rate = migration_fraction / migration_time # miogration rate for ms given the fraction and length
migration_start = tea_time / (4 * Ne) # time at which migrattion star given the donor et the length of the migration
migration_end = migration_start + migration_time # time at which the migration stop
mutation_rate = len_locus * 4 * Ne * mu

# R source code variables
Head = "library('ape') \n" + "library('coala') \n" + "library('phyclust') \n" + "library('phangorn') \n" + "activate_ms(priority = 600) \n\n"
Pop = "pop <- c(%s) \n\n" % ', '.join(str(x) for x in ext_lineages)
D_R = "name_donor <- '%s' \nname_recip <- '%s' \n\n" % (pop_donor.name, pop_recip.name)
Coal = "model <- coal_model(sample_size = c(%s), loci_number = %s, loci_length = 1, ploidy = %s) + \n" % (', '.join(str(x) for x in ext_lineages), len_locus, ploidy)
# Mutation = "feat_mutation(rate = %s, model = 'IFS', fixed_number = FALSE, locus_group = 'all') + \n" % mutation_rate
Mutation = "feat_mutation(1, fixed_number = TRUE, locus_group = 'all') + \n"
Migration_starts = "feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    migration_rate, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_start)
Migration_ends = "feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_end)
Stat_sum = "sumstat_trees() \n"
True_migration = "\n# donor: " + pop_donor.name + " -- > recipient: " + pop_recip.name + ' at ' + str(n_generation_to_root - tea_time) + "\n"

# R source command output
# with open(os.path.join(args.output,'ms_command.R'), "w") as output:
with open('ms_command_ghost.R', "w") as output:
    output.write("".join([Head, Pop, D_R, Coal, Mutation, Merge, Migration_starts, Migration_ends, Stat_sum, True_migration]))

ta = t.get_common_ancestor(trio)
tp = ta.copy()

n_generation_to_root = int(tp.name.split("@")[1])


Head = D_R = Pop = Coal = Merge = Migration_starts = Migration_ends = Stat_sum = True_migration = str() # coal_model extended parameters
count = 1
for i in tp.traverse('postorder'):
    old=i.name
    if i.is_leaf():
        ext_lineages.append(1)
        i.name = "ali>" + str(count) + "@0"
        count += 1
    else:
        pop_d = str(gene_n(i.get_descendants()[0])).split("_")[0]
        pop_r = str(gene_n(i.get_descendants()[1])).split("_")[0]
        pop_g = int(n_generation_to_root - i.get_distance(tp))
        i.name = pop_d + "_" + pop_r + "@" + str(pop_g)
        pop_t =  pop_g / (4 * Ne)
        Merge += "feat_pop_merge(%s, %s, %s) + " % (str(pop_t), pop_r, pop_d) + "\n"


pop_recip = tp&"ali>1@0"
pop_donor = tp&"ali>2@0"

tea_time = int(
    (max(int(pop_recip.name.split("@")[1]), int(pop_donor.name.split("@")[1]))+
    min(int(pop_recip.up.name.split("@")[1]), int(pop_donor.up.name.split("@")[1])))/2)


migration_time = 1 / (4 * Ne) # time of the migration in 4Ne
migration_rate = migration_fraction / migration_time # miogration rate for ms given the fraction and length
migration_start = tea_time / (4 * Ne) # time at which migrattion star given the donor et the length of the migration
migration_end = migration_start + migration_time # time at which the migration stop
mutation_rate = len_locus * 4 * Ne * mu



Head = "library('ape') \n" + "library('coala') \n" + "library('phyclust') \n" + "library('phangorn') \n" + "activate_ms(priority = 600) \n\n"

D_R = "name_donor <- '%s' \nname_recip <- '%s' \n\n" % (pop_donor.name, pop_recip.name)

Coal = "model <- coal_model(sample_size = c(1, 1, 1), loci_number = %s, loci_length = 1, ploidy = %s) + \n" % (len_locus, ploidy)

Mutation = "feat_mutation(1, fixed_number = TRUE, locus_group = 'all') + \n"

Migration_starts = "feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    migration_rate, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_start)

Migration_ends = "feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_end)

Stat_sum = "sumstat_trees() \n"

True_migration = "\n# donor: " + pop_donor.name + " -- > recipient: " + pop_recip.name + ' at ' + str(n_generation_to_root - tea_time) + "\n"


with open('ms_command_sampled.R', "w") as output:
    output.write("".join([Head, Pop, D_R, Coal, Mutation, Merge, Migration_starts, Migration_ends, Stat_sum, True_migration]))



#GNU Terry Pratchett
