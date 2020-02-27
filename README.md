# abba_baba


Running a simulation:

    run.sh  /PATH/Singularity.sif Simulation_output seed Parameters/zombi_parameters Parameters/sim_parameters


Exemple:

    run.sh  sim.sif sim_1 8059 Parameters/zombi_parameters Parameters/sim_parameters


If seed = 0, then a random seed is used for each set, else the same seed is used for the tree simulation and the rest of the pipeline.
