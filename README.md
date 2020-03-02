# abba_baba


Running a simulation:

    run.sh  /PATH/Singularity.sif Simulation_output seed Parameters/zombi_parameters Parameters/sim_parameters


Exemple:

    run.sh  sim.sif sim_1 8059 Parameters/zombi_parameters Parameters/sim_parameters


If seed = 0, then a random seed is used for each set, else the same seed is used for the tree simulation and the rest of the pipeline.


You can use singularirty to build a container with all the required dependency for this work.
You need to be sudo to perform this task.

    sudo singularity build singularity.sif ~/GitHub/abba_baba/singularity/singularity.def


ZOMBi version in the singularity container is up to the commit 21d562cc380fafca58a7068913f4a17f52c60451
If you are using an up to date ZOMBI you will have to use ZOMBI/Partameters/SpeciesTreeParameters.tsv instead of abba_baba/Parameters/zombi_parameters.
This can lead to unforeseen bugs!


<!--
for i in `seq 2999 1 3099`; do
 sed "s/aaaa/$i/g" run_slurm.sh > temp
 sbatch temp
 rm temp
done
 -->
