# Theo

#run ms simulation for ABBA BABA


for i in `seq 1 1 9999999`; do
  if [ ! -d sim_$i]; then
    dir=sim_$i
    python ZOMBI/ZOMBI.py T ZOMBI/Parameters/SpeciesTreeParameters.tsv $dir
    python python3/build_ms_command.py $dir/T/CompleteTree.nwk -p ~/GitHub/abba_baba/shell/Parameters -o $dir -v
    # seaview $dir/spe_tree &
    Rscript R/ms_simulation.R $dir/ms_command.R
    break
    stop
  fi
done



for i in `seq 1 1 99999999`; do
  dir=sim_$i
  rm -rf $dir
  zombi T SpeciesTreeParameters.tsv $dir
  # rm -rf $dir/SAM*
  # py ~/GitHub/ZOMBI/SpeciesSampler.py n 20 $dir  # seaview $dir/T/CompleteTree.nwk &

  py ~/GitHub/abba_baba/python3/build_ms_command.py $dir/T/CompleteTree.nwk -p ~/GitHub/abba_baba/shell/Parameters -o $dir -v
  # seaview $dir/spe_tree &
  Rscript ~/GitHub/abba_baba/R/ms_simulation.R $dir/ms_command.R
done


##
# singularity

# building from def
sudo singularity build sim.sif singularity.def

singularity shell singularity.sif --cleanenv --contain --no-home

ssh -AX ubuntu@134.214.213.82

rsync -a singularity.sif tricou@pbil-deb.univ-lyon1.fr:/beegfs/home/tricou/singularity_image

rsync -a singularity.sif ubuntu@134.214.213.82:


########################

dir=sim_1
cp abba_baba/Parameters $dir/Param_$dir
python /usr/local/ZOMBI/ZOMBI.py T $dir/Param_$dir $dir
python abba_baba/python3/build_ms_command.py $dir/T/CompleteTree.nwk -p $dir/Param_$dir -o $dir -v
Rscript abba_babaR/ms_simulation.R $dir/ms_command.R
















dir=sim_1
mkdir sim_1
cp shell/Parameters $dir/Param_$dir

zombi T zombi_parameters $dir

python python3/build_ms_command.py $dir/T/CompleteTree.nwk -p shell/Parameters -o $dir -v
# seaview $dir/spe_tree &
Rscript R/ms_simulation.R $dir/ms_command.R shell/Parameters






singularity run singularity.sif --R R/ms_simulation.R
