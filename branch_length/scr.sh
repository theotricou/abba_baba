#!/bin/bash
#Theo

# branch length paper

# mosquito type data


"""
simulate a species tree with ghost lineges (unsam + ext)

selecte a group of 3 lineages with (A,B)C)
and ancestor(A,B,C) coalescent aproximatively at 10 % of the size spe tree

simulate 1 (one) introgression from ghost to A in Complete spe tree

simulate 1 (one) introgression from ghost to A in Sampled spe tree

Compare branche length of all 3 lineages topology
"""

rm -rf test tree && py ~/GitHub/Zombi/Zombi.py T SpeciesTreeParameters.tsv test
# sss test/T/CompleteTree.nwk &
py ~/GitHub/unknown_diversity/hgt_ghost/ladderize_tree.py test/T/CompleteTree.nwk

# sss tree &
sss test/T/CompleteTree.nwk &




# GNU Pratchett
