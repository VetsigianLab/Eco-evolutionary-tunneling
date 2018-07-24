1) To generate Fig. 1C

   >> load simulation_output_Fig5c %this the recorded trajectory
   >> plot_spatial_sim %this plots a simulated trajectory

2) To simulate a trajectory see set_up_and_run_spatial_simulation.m

   Basically, set up a simulation structure and run

   evolve_RS_MutOne_Space('file_name', simulation)


3) ESS_3_153_mu4_r10.mat contains a non-spatial simulation that leads to an ESS. 
   Its final state and parameters are used for teh spatial simulation. 