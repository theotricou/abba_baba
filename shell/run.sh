# Theo

# run ms simulation for ABBA BABA



cp ~/GitHub/ZOMBI/Parameters/SpeciesTreeParameters.tsv .

zombi T SpeciesTreeParameters.tsv tree_simulation

seaview tree_simulation/T/CompleteTree.nwk

cp tree_simulation/T/CompleteTree.nwk tree

py ~/GitHub/abba_baba/python3/build_ms_command.py tree_simulation/T/CompleteTree.nwk ~/GitHub/abba_baba/shell/Parameters

Rscript ~/GitHub/abba_baba/R/ms_simulation.R ms_command.R Parameters
