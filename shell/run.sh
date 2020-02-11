# Theo

# run ms simulation for ABBA BABA


path=$PWD

cp ~/GitHub/ZOMBI/Parameters/* .
zombi T SpeciesTreeParameters.tsv tree_sim

pyton3 build_ms_command.py tree_sim/T/CompleteTree.nwk

Rscript ms_simulation.R ms_command.txt
