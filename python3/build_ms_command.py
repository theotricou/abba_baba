#! /usr/bin/env python
# Theo

# build ms command file from a phylogeny tree

import os
import sys
from ete3 import Tree as tr
import random
import argparse

# exp_n = 0
# dir = "Sim_" + str(exp_n)
# while os.path.isdir(dir) == True:
#     dir = "Sim_" + str(exp_n)
#     exp_n += 1

dir = "Simulation"
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('tree', help ='A required phylogenetic tree argument')
parser.add_argument('-p', '--parameters',action="store_true", help ='An optionnal parameters file argument')
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("output", type=str, nargs="?", default= dir, help="Name of the simulation folder")
args = parser.parse_args()

if os.path.isdir(args.output):
    sys.exit("Output folder already present in experiment folder. Please, remove previous existing data to proceed.")
else:
    os.mkdir(args.output)


def read_param(string):
        with open(parameters_file) as f:
            for line in f:
                if line[0] == "#" or line == "\n":
                    continue
                if string + "=" in line:
                    parameter, value = line.strip().split("=")
                    break
        return(float(value))

def any_descendant_alive(node):
    # test if the recipient of a gene flow has any descendant, if not no gene flow can be detecte
    if node == False:
        return(False)
    else:
        for i in node.get_leaves():
            if int(ext_lineages[int(i.name) - 1]) == 1:
                return(True)
                break
            else:
                return(False)

def search_alive_togethere(tree, node):
    "Finds nodes alive at the same time as the donor"
    matches = [node.get_sisters()[0]]
    for n in tree.iter_descendants('postorder'):
        if n == node or n == node.get_sisters()[0]:
            pass
        else:
            if node.gen <= n.gen < node.up.gen or node.gen < n.up.gen <= node.up.gen:
                matches.append(n)
    return(matches)

# general ms parameters
# if the parameters file existe the following default parameters are overruled

if args.parameters == False:
    Ne = 100000
    mu = 1e-7
    len_locus = 1
    ploidy = 1
    n_genaration_to_root = 1600000
else:
    parameters_file = args.parameters
    Ne = read_param("NE")
    mu = read_param("MU")
    len_locus = read_param("LOCI_LENGTH")
    ploidy = read_param("PLOIDY")
    n_genaration_to_root = read_param("N_GENERATION")


t = tr(args.tree, format = 1) # read phylo tree

# variables for the outputing of the source code for R
ext_lineages = [] # list used for coal_model population
Head = D_R = Pop = Coal = Merge = Migration_starts = Migration_ends = Stat_sum = str() # coal_model extended parameters

# build the merge part of the ms command
# an ext_lineages variable is crated to nkow which population ms samples
t.gen = n_genaration_to_root
name = 1
for i in t.traverse('postorder'):
    i.add_features(gen = 0)
    if i.is_leaf():
        i.name = name
        name += 1
        if i.get_distance(t) < t.get_farthest_leaf()[1] * 0.99:
            ext_lineages.append("0")
            i.gen = round(((t.get_farthest_leaf()[1] - i.get_distance(t)) * n_genaration_to_root / t.get_farthest_leaf()[1]), -3) / (4 * Ne)
        else:
            ext_lineages.append("1")
    else:
        i.gen = round(((t.get_farthest_leaf()[1] - i.get_distance(t)) * n_genaration_to_root / t.get_farthest_leaf()[1]), -1) / (4 * Ne)
        i.name = str(i.get_descendants()[0].name).split("_")[0] + "_" + str(i.get_descendants()[1].name).split("_")[0] # name used in feat_pop_merge
        Merge += "feat_pop_merge(%s, %s, %s) + " % (i.gen, str(i.get_descendants()[1].name).split("_")[0], str(i.get_descendants()[0].name).split("_")[0]) + "\n"


# migration Parameters
migration_generation = 10 # number of generation during which the population migrate
migration_fraction = 0.1 # fraction of the donor population that migrate in recipient population, during migration_generation generations
migration_time = migration_generation / (4 * Ne) # time of the migration in 4Ne
migration_rate = migration_fraction / migration_time # miogration rate for ms given the fraction and length


range_migration = [0,0]
# randomly choose a donor lineage and a recipient (this recipient need to have a least one extant descendant)
pop_donor = pop_recip = False
while any_descendant_alive(pop_recip) == False or migration_generation > ((range_migration[0] * 4 * Ne) - (range_migration[1] * 4 * Ne) + (migration_generation / 2)):
    pop_donor = random.choice(t.get_descendants())
    pop_recip = random.choice(search_alive_togethere(t, pop_donor))
    range_migration = [min(pop_donor.get_ancestors()[0].gen, pop_recip.get_ancestors()[0].gen), max(pop_donor.gen, pop_recip.gen)] # range during which the migration is possible, migration can't past a population merge


migration_start = round(min(range_migration), 9) + ((migration_generation / 2) / (4 * Ne)) # time at which migrattion star given the donor et the length of the migration
migration_end = round(migration_start + migration_time, 9) # time at which the migration stop
mutation_rate = len_locus * 4 * Ne * mu

# R source code variables
Head = "library('ape') \n" + "library('coala') \n" + "library('phyclust') \n" + "library('phangorn') \n " + "library('data.table') \n " + "activate_ms(priority = 600) \n"


Pop = "pop <- c(%s) \n" % ', '.join(ext_lineages)

D_R = "name_donor <- '%s' \n name_recip <- '%s' \n" % (pop_donor.name, pop_recip.name)

Coal = "model <- coal_model(sample_size = c(%s), loci_number = %s, loci_length = 1, ploidy = %s) + \n" % (', '.join(ext_lineages), len_locus, ploidy)

Mutation = "feat_mutation(rate = %s, model = 'IFS', fixed_number = FALSE, locus_group = 'all') + \n" % mutation_rate

Migration_starts = "feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (migration_rate, str(pop_donor.name).split("_")[0], str(pop_recip.name).split("_")[0], migration_start)
Migration_ends = "feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (str(pop_donor.name).split("_")[0], str(pop_recip.name).split("_")[0], migration_end)

Stat_sum = "sumstat_seg_sites() + sumstat_trees() \n"

# tree output, the basic species tree is output in spe_tree
# in extended_spe_tree the coresponding generation separating each node is output in an extended newik format
t.write(outfile = os.path.join(args.output,'extended_spe_tree'), format=1, features=["gen"], format_root_node=True)
t.write(outfile = os.path.join(args.output,'spe_tree'), format=1, format_root_node=True)

# R source command output
with open(os.path.join(args.output,'ms_command.R'), "w+") as output:
    output.write("".join([Head, Pop, D_R, Coal, Mutation, Merge, Migration_starts, Migration_ends, Stat_sum]))

if args.verbose:
    print(pop_donor.name, pop_recip.name)
