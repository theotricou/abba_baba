# Theo

#run ms simulation for ABBA BABA

rm summary_simulation

cp ~/GitHub/zombi/Parameters/SpeciesTreeParameters.tsv


for i in `seq 2 1 5`; do
  dir=sim_$i
  rm -rf $dir
  zombi T SpeciesTreeParameters.tsv $dir
  rm -rf $dir/SAM*
  py ~/GitHub/ZOMBI/SpeciesSampler.py n 20 $dir  # seaview $dir/T/CompleteTree.nwk &

  py ~/GitHub/abba_baba/python3/build_ms_command.py $dir/T/CompleteTree.nwk -p ~/GitHub/abba_baba/shell/Parameters -o $dir -v
  # seaview $dir/spe_tree &
  Rscript ~/GitHub/abba_baba/R/ms_simulation.R $dir/ms_command.R
done


##
# singularity

# create sandbox dir named ubuntu
sudo singularity build --sandbox ubuntu/ library://ubuntu


# building from def
sudo singularity build singularity.sif test_singularity.def
