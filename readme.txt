The *.0  files are input Spparks files.

Execute the follwing command to run the Spparks input files (*.0)
	$ mpirun -np 1 spk_redsky < *.0

spk_redsky is an executable file generated while building Spparks. The instruction for building spparks executable can be found at the website,  https://spparks.github.io/doc/Section_start.html


The source codes for Spparks to model NK signaling for different models and NK cell signaling are attached in the respective folders.

	1. model1_act_sig_spparks_src_codes (source codes to stimulate Model 1 in presence of activating ligands; the rates are estimated by PSO)

	2. model2_act_sig_spparks_src_codes (source codes to stimulate Model 2 in presence of activating ligands; the rates are estimated by PSO)

	3. model3_act_sig_spparks_src_codes (source codes to stimulate Model 3 in presence of activating ligands; the rates are estimated by PSO)

	4. model2_inb_sig_spparks_src_codes (source codes to stimulate Model 2 in presence of both activating and inhibitory ligands (inhibitory ligands are distributed heterogeneously e.g.ring); the rates are estimated by PSO)

	5. model3_inb_sig_spparks_src_codes (source codes to stimulate Model 3 in presence of both activating and inhibitory ligands (inhibitory ligands are distributed heterogeneously e.g.ring); the rates are estimated by PSO)




model2_PSO (In-silco framework to estimate the model parameters via Particle Swarm Optimization, PSO; model under consideration is Model 2 )

The following files are included in the above folder:

	1. run_pso.py (Code to calculate the cost function with estimated parameters on logarithmic scale written in an output file)

	2. setup_full_run.py (Code to create Spparks input file for a given set of parameters estimated via PSO; model under consideration is Model 2)

	3. voxel_size.py (Code to calculate the number of NKG2D within the simulation box of size 15 x 15 micron for a given set of parameters estimated via PSO)

	4. std_density.py (Code to calculate the mean and standard deviation of NKG2D molecules for a given set of parameters estimated via PSO)

	5. pair_corr.py (Code to calculate the two point correlation function for NKG2D molecules within the simulation box of size 15 x 15 micron for a given set of parameters estimated via PSO)
 
	7. S4_std_density.py (Code to calculate the mean and standard deviation of coarse-grained data for Dap10-mCherry from the TIRF images from Ref. [1] at time t = 60 sec)

	8. thres_S4_2_avg_voxel.py (Code to calculate the two point correlation function for coarse-grained data for Dap10-mCherry from the TIRF images from Ref. [1] at time t = 60 sec)
	
	9. Pyswarm/pso.py (code for asynchronous Particle Swarm Optimization (PSO) with constraints - modified from the standard PSO algorithm available at the link https://github.com/tisimst/pyswarm )

	10. sites.org_S4_1_60.60 (Input file: extracted intensity data for Dap10-mCherry corresponding to the region of interest from the TIRF images from Ref. [1] at time t = 60 sec )

	

Ref[1]: Thushara P. Abeyweera, Ernesto Merino, Morgan Huse. Inhibitory signaling blocks activating receptor clustering and induces cytoskeletal retraction in natural killer cells. J Cell Biol, 192 (4): 675–690, (2011).