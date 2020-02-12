#! /usr/bin/env python
# Theo

# build ms command file from a phylogeny tree

import os
import sys
from ete3 import Tree as tr
import random

# tree = sys.argv[1]
tree = "tree"

t = tr(tree, format = 1)

name = 1
Ne = 100000
mu = 1e-7
n_loci = 1
len_locus = 1
ploidy = 1
n_genaration_to_root = 1600000
t.gen = n_genaration_to_root

ext_lineages = [] # list used for coal_model
Head = D_R = Pop = Coal = Merge = Migration_starts = Migration_ends = Stat_sum = str()

for i in t.traverse('postorder'):
    i.add_features(gen = 0)
    if i.is_leaf():
        i.name = name
        name += 1
        if i.get_distance(t) < t.get_farthest_leaf()[1] * 0.99:
            ext_lineages.append("0")
        else:
            ext_lineages.append("1")
    else:
        i.gen = round(((t.get_farthest_leaf()[1] - i.get_distance(t)) * n_genaration_to_root / t.get_farthest_leaf()[1]), -3) / (4 * Ne)
        i.name = str(i.get_descendants()[0].name).split("_")[0] + "_" + str(i.get_descendants()[1].name).split("_")[0] # name used in feat_pop_merge
        Merge += "feat_pop_merge(%s, %s, %s) + " % (i.gen, str(i.get_descendants()[1].name).split("_")[0], str(i.get_descendants()[0].name).split("_")[0]) + "\n"

# randomly choose a donor lineage and a recipient (this recipient need to have a least one extant descendant)
def any_descendant_alive(node):
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
    matches = []
    for n in tree.iter_descendants('postorder'):
        if n == node:
            pass
        else:
            if node.gen <= n.gen < node.up.gen or node.gen < n.up.gen <= node.up.gen:
                matches.append(n)
    return(matches)


migration_generation = 10000 # number of generation during which th population migrate
migration_fraction = 0.1 # fraction of the donor population that migrate in recipient population during migration_generation generation
migration_time = migration_generation / (4 * Ne) # time of the migration in 4Ne
migration_rate = migration_fraction / migration_time # miogration rate for ms given the fraction and length

pop_donor = pop_recip = False
#pop_donor = random.choice(t.get_descendants())
while any_descendant_alive(pop_recip) == False or migration_generation > ((range_migration[0] * 4 * Ne) - (range_migration[1] * 4 * Ne) + (migration_generation / 2)):
    pop_donor = random.choice(t.get_descendants())
    pop_recip = random.choice(search_alive_togethere(t, pop_donor))
    range_migration = [min(pop_donor.get_ancestors()[0].gen, pop_recip.get_ancestors()[0].gen), max(pop_donor.gen, pop_recip.gen)] # range during which the migration is possible, migration can't past a population merge


migration_start = round(min(range_migration), 9) + ((migration_generation / 2) / (4 * Ne)) # time at which migrattion star given the donor et the length of the migration
migration_end = round(migration_start + migration_time, 9) # time at which the migration stop
mutation_rate = len_locus * 4 * Ne * mu


Head = "require('coala') \n" + "require('phyclust') \n" + "activate_ms(priority = 600) \n"

Pop = "pop <- c(%s) \n" % ', '.join(ext_lineages)

D_R = "name_donor <- '%s' \n name_recip <- '%s' \n" % (pop_donor.name, pop_recip.name)

Coal = "model <- coal_model(sample_size = c(%s), loci_number = %s, loci_length = %s, ploidy = %s) + \n" % (', '.join(ext_lineages), n_loci, len_locus, ploidy)

Mutation = "feat_mutation(rate = %s, model = 'IFS', fixed_number = FALSE, locus_group = 'all') + \n" % mutation_rate


Migration_starts = "feat_migration(%s, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (migration_rate, str(pop_donor.name).split("_")[-1], str(pop_recip.name).split("_")[-1], migration_start)
Migration_ends = "feat_migration(0, pop_from = %s, pop_to = %s, symmetric = FALSE, time = %s, locus_group = 'all') + \n" % (str(pop_donor.name).split("_")[0], str(pop_recip.name).split("_")[0], migration_end)

Stat_sum = "sumstat_seg_sites() + sumstat_trees() \n"

t.write(outfile = 'Species_tree', format=1, features=["gen"], format_root_node=True)

with open('ms_command.R', "w+") as output:
    output.write("".join([Head, Pop, D_R, Coal, Mutation, Merge, Migration_starts, Migration_ends, Stat_sum]))

print(pop_donor.name, pop_recip.name)
