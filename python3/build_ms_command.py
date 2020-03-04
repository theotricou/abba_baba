#! /usr/bin/env python
# Theo

# build ms command file from a phylogeny tree

import os
import sys
from ete3 import Tree as tr
import random
import argparse
import numpy as np

print("\nBuilding ms command")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('tree', type = str, help = 'A phylogenetic tree')
parser.add_argument('-o', '--output', type = str, nargs = "?", default = "Simulation", help = "Name of the simulation folder. Default is 'Simulation'.")
parser.add_argument('-p', '--parameters', type = str, nargs = "?", default = False, help ='An optionnal parameters file argument')
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
args = parser.parse_args()

def gene_t(node):
    return(str(node.name).split('@')[1])

def gene_n(node):
    return(str(node.name).split('@')[0])

def read_param(string):
        with open(parameters_file) as f:
            for line in f:
                if line[0] == "#" or line == "\n":
                    continue
                if string + "=" in line:
                    parameter, value = line.strip().split("=")
                    break
        return(float(value))


def tree_new_dist(tree, n_genaration_to_root):
    multi = n_genaration_to_root / tree.get_farthest_leaf()[1]
    for i in tree.traverse():
        if i.is_leaf():i.dist = int(round(i.dist * multi))
        else : i.dist = int(round(i.dist * multi) + 1)
    for i in tree:
        if i.is_leaf():
            if i.get_distance(tree) > n_genaration_to_root:
                i.dist = int(i.dist - (i.get_distance(tree) - n_genaration_to_root))
    return(tree)

def alive_at_time(tree, time):
    results = []
    for i in tree.iter_descendants():
        if i.up.dist < time & time < i.get_distance(t):
            results.append(i)
    return(results)

def any_descendant_alive(tree, node):
    # test if the recipient of a gene flow has any descendant, if not no gene flow can be detecte
    if node == False:
        return(False)
    else:
        if "@0" in node.name:
            return(True)
        else:
            for i in node.get_leaves():
                if "@0" in i.name:
                    return(True)
                    break
                else:
                    return(False)

# general ms parameters,if the parameters file existe the following default parameters are overruled
if args.parameters == False:
    Ne = 100000
    mu = 1e-7
    len_locus = 1
    ploidy = 1
    n_genaration_to_root = 1000000
else:
    parameters_file = args.parameters
    Ne = read_param("NE")
    mu = read_param("MU")
    len_locus = read_param("LOCI_LENGTH")
    ploidy = read_param("PLOIDY")
    n_genaration_to_root = read_param("N_GENERATION")
    # os.system('cp %s %s' % (args.parameters, args.output))
    if not read_param("SEED") == 0:
        random.seed(int(read_param("SEED")))
        np.random.seed(int(read_param("SEED")))


t = tr(args.tree, format = 1) # read phylo tree
t.dist = 0
t = tree_new_dist(t, n_genaration_to_root)
sample = []
if os.path.isfile(os.path.join(*args.tree.split('/')[0:-2], "SAMPLE_1/SampledSpeciesTree.nwk")):
    ts = tr(os.path.join(*args.tree.split('/')[0:-2], "SAMPLE_1/SampledSpeciesTree.nwk"), format = 1) # read phylo extant tree
    for i in ts: sample.append(i.name)
else:
    te = tr(os.path.join(*args.tree.split('/')[0:-1], "ExtantTree.nwk"), format = 1) # read phylo extant tree
    for i in te: sample.append(i.name)

# variables for the outputing of the source code for R
Head = D_R = Pop = Coal = Merge = Migration_starts = Migration_ends = Stat_sum = True_migration = str() # coal_model extended parameters
ext_lineages = []
count = 1
for i in t.traverse('postorder'):
    if i.is_leaf():
        if i.name in sample:
            ext_lineages.append(1)
        else:
            ext_lineages.append(0)
        temp_ = str(int(n_genaration_to_root - i.get_distance(t)))
        i.name = str(count) + "@" + temp_
        count += 1
    else:
        pop_d = str(gene_n(i.get_descendants()[0])).split("_")[0]
        pop_r = str(gene_n(i.get_descendants()[1])).split("_")[0]
        pop_g = int(n_genaration_to_root - i.get_distance(t))
        i.name = pop_d + "_" + pop_r + "@" + str(pop_g)
        pop_t =  pop_g / (4 * Ne)
        Merge += "feat_pop_merge(%s, %s, %s) + " % (str(pop_t), pop_r, pop_d) + "\n"

# randomly choose a donor lineage and a recipient (this recipient need to have a least one extant descendant)
if args.verbose:
    print("\nPicking migration duo. Possible bottleneck!")

# all_node = t.get_descendants()
pop_donor = pop_recip = False
while any_descendant_alive(t, pop_recip) == False:
    time_mig = int(np.random.uniform(0, n_genaration_to_root))
    all_node =  alive_at_time(t, time_mig)
    pop_donor, pop_recip = np.random.choice(all_node, size = 2, replace = False)

# migration Parameters
migration_generation = 1 # number of generation during which the population migrate
migration_fraction = 0.1 # fraction of the donor population that migrate in recipient population, during migration_generation generations
migration_time = migration_generation / (4 * Ne) # time of the migration in 4Ne
migration_rate = migration_fraction / migration_time # miogration rate for ms given the fraction and length
migration_start = time_mig / (4 * Ne) # time at which migrattion star given the donor et the length of the migration
migration_end = migration_start + migration_time # time at which the migration stop
# mutation parameter
mutation_rate = len_locus * 4 * Ne * mu

# R source code variables
Head = "library('ape') \n" + "library('coala') \n" + "library('phyclust') \n" + "library('phangorn') \n" + "activate_ms(priority = 600) \n\n"
Pop = "pop <- c(%s) \n\n" % ', '.join(str(x) for x in ext_lineages)
D_R = "name_donor <- '%s' \nname_recip <- '%s' \n\n" % (pop_donor.name, pop_recip.name)
Coal = "model <- coal_model(sample_size = c(%s), loci_number = %s, loci_length = 1, ploidy = %s) + \n" % (', '.join(str(x) for x in ext_lineages), len_locus, ploidy)
Mutation = "feat_mutation(rate = %s, model = 'IFS', fixed_number = FALSE, locus_group = 'all') + \n" % mutation_rate
Migration_starts = "feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    migration_rate, str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_start)
Migration_ends = "feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (
    str(gene_n(pop_donor)).split("_")[0], str(gene_n(pop_recip)).split("_")[0], migration_end)
Stat_sum = "sumstat_seg_sites() + sumstat_trees() \n"
True_migration = "\n# donor: " + pop_donor.name + " -- > recipient: " + pop_recip.name + "\n"

t.write(outfile = os.path.join(args.output,'spe_tree'), format=1, format_root_node=True)

# R source command output
with open(os.path.join(args.output,'ms_command.R'), "w") as output:
    output.write("".join([Head, Pop, D_R, Coal, Mutation, Merge, Migration_starts, Migration_ends, Stat_sum, True_migration]))

if args.verbose:
    print(True_migration)
