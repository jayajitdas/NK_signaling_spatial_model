------------------
Folder: spparks_src_codes:
------------------

Each sub-directories contains source codes ( *.cpp and *.h files for SPPARKS) for different models for NK cell signaling.

	1. model1_act_spparks_src_codes (codes to Model 1 for simulation by activating ligands; the rates are estimated by PSO)
	2. model2_act_spparks_src_codes (codes to Model 2 for simulation by activating ligands; the rates are estimated by PSO)
	3. model3_act_spparks_src_codes (codes to Model 3 for simulation by activating ligands; the rates are estimated by PSO)
	4. model2_inb_spparks_src_codes (codes to Model 2 for simulation by both activating and inhibitory ligands; 
	   the rates are estimated by PSO)
	5. model3_inb_spparks_src_codes (codes to Model 2 for simulation by both activating and inhibitory ligands; 
	   the rates are estimated by PSO)
	6. model1_act_no_clus_spparks_src_codes (codes to Model 1 for simulation by activating ligands without NKG2D cluster formation; 
	   the rates are estimated by PSO)
	7. model2_act_no_clus_spparks_src_codes (codes to Model 2 for simulation by activating ligands without NKG2D cluster formation; 
	   the rates are estimated by PSO)

Step1: Download Spparks from https://sjplimp.github.io//download.html
----------------------------------------------------------------------

To run a specific model, add the files (*.cpp and *.h) from the corresponding sub-directory to the "src" folder in the Spparks downloaded folder.

For eg: If you want to simulate  Model 2 in presence of activating ligands only. Copy all *.cpp and *.h files from model2_act_spparks_src_codes 
to the src folder of the SPPARKS downloaded folder.

Step2: Build sparks executable: 
--------------------------------
# go to directiry ~spparks_downloaded_folder/src
# run following command:
  make redsky
# an executable file with the name spk_redsky will be generated.




-------
Models:
-------

Each of these sub-directories contains a model you can run using SPPARKS. Each model has an input script (in.*), input read_file (sites.30.30), 
a bash script (run.sh) and produces a log file (log.*) and output files (sites.*.*) when they run. 

The sites.*.* files produced by the model runs can be converted to *.vtk to visualize using ParaView-5.7.0

Step3: How to run and visualize a particular model (For eg: Model 2 in presence of activating ligands):
-------------------------------------------------------------------------------------------------------
cd model2_act	# go to the work directory
cp ../~spparks_downloaded_folder/src/spk_redsky . #copy the executable file spk_redsky (from Step2) 
						  build using ???model2_act_sig_spparks_src_codes??? from
						 ~spparks_downloaded_folder/src folder to the current
						  directory
sbatch run.sh	# To run the code

More details about Spparks can be found at the website,  https://spparks.github.io/doc/Section_start.html

These are the models in the various sub-directories:

	1. model1_act (Model 1 in presence of activating ligands; 
	   the rates are estimated by PSO)
	2. model1_act_no_clus (Model 1 in presence of activating ligands without clustering of NKG2D; 
           the rates are estimated by PSO)
	3. model2_act (Model 2 in presence of activating ligands; 
	   the rates are estimated by PSO)
	4. model2_act_no_clus (Model 2 in presence of activating ligands without clustering of NKG2D; 
	   the rates are estimated by PSO)
	5. model2_inb (Model 2 in presence of both activating and inhibitory ligands 
	   (inhibitory ligands are distributed heterogeneously e.g.ring); 
	    the rates are estimated by PSO)
	6. model3_act (Model 3 in presence of activating ligands ; the rates are estimated by PSO)
	7. model3_inb(Model 3 in presence of both activating and inhibitory ligands 
	   (inhibitory ligands are distributed heterogeneously e.g.ring); 
	   the rates are estimated by PSO)

Step4: Convert sites.*.* output files from each of the sub-directories to *.vtk for visualization using ParaView
-------------------------------------------------------------------------------------------------------
Run the following command:

for i in `seq 0 180`;do if [ -f sites.0.$i ];then echo "sites.0.$i => ar_only.$i.vtk"; if [ ! -f ar_only.$i.vtk ]; then ../sites2vtk_all_ar.py <sites.0.$i > ar_only.$i.vtk; fi; fi;done


The above command will convert each sites.*.* output fils from the model runs to corresponding ar_only.*.vtk.
The resulting ar_only.*.vtk files can be used to visualize the model data using Paraview-5.7.0 (https://www.paraview.org/download/)





---------------------
Parameter_estimation :
---------------------

In-silco framework to estimate the model parameters via Particle Swarm Optimization (PSO):

Each sub-directory contains files to run PSO in parallel to estimate the parameters for a model under consideration (Model_1 and Model_2).  

	1. run.sh (bash script to run the In-silico framework)
	2. run_pso.py (calculate the cost function with estimated parameters on logarithmic scale)
	3. setup_full_run.py (creates Spparks input file for the Model given for a set of parameters
	   estimated via PSO; model under consideration is Model 1 or Model 2)
	4. voxel_size.py (calculate the number of NKG2D within the simulation box of size 15 x 15 micron 
	   for a given set of parameters estimated via PSO)
	5. std_density.py (calculate the mean and standard deviation of NKG2D molecules for a given set 
	   of parameters estimated via PSO)
	6. pair_corr.py (calculate two point correlation function for NKG2D molecules within the 
	   simulation box of size 15 x 15 micron for a given set of parameters estimated via PSO)
	7. S4_adjusted_pair_corr.py (calculate two point correlation function of coarse-grained data for 
	   Dap10-mCherry from the TIRF images S4 from Ref. [1] at time t = 60 sec)
	8. S4_std_density.py (calculate mean and standard deviation of coarse-grained data for 
	   Dap10-mCherry from the TIRF images S4 from Ref. [1] at time t = 60 sec)
	9. thres_S4_2_avg_voxel.py (calculate two point correlation function for coarse-grained data 
	   for Dap10-mCherry from the TIRF images S4 from Ref. [1] at time t = 60 sec)	
	10. Pyswarm/pso.py (code for asynchronous Particle Swarm Optimization (PSO) with constraints 
	   - modified from the standard PSO algorithm available at the link 
	   https://github.com/tisimst/pyswarm )
	11. sites.org_S4_1_60.60 (Input file: extracted intensity data for Dap10-mCherry 
	    corresponding to the region of interest from the TIRF images S4 from Ref. [1] 
	    at time t = 60 sec )
	12. spk_redsky (spparks executable file build for the respective Model)
	13. spk_redsky_no_clus (spparks executable file build for case of Model 2 without NKg2D clustering)

To get output run the following command from a bash script:
-----------------------------------------------------------
cd parameter_estimation
cd Model_1 # For eg: run parameter estimation fro Model 1
sbatch run.sh # to run the framework


This will generate a series of folders with prefix ???md5_*??? and a output file (slurm.*.out). Each row of the output file (slurm.*.out) contains:
 	1. name of the folder (md5_*)containing Spparks input and output files for a given set of parameters.
	2. list of respective set of parameters estimated by PSO.  
	3. cost function value corresponding to the respective set of parameters estimated by PSO. 

Note:
-----
$ Input file ???sites.org_S4_1_60.60??? can be replaced with the intensity profile of the image one wants to estimate the parameters for using our Model 2.

$To change the model, replace the both setup_full_run.py and spk_redsky with executable file build for the respective model.





---------------------------------
PSO uncertainties quantification:
---------------------------------


To quantify the uncertainties in the estimated parameters we used Density Peak Clustering algorithm in Ref.[2]. The folder contains to files:

	1. norm_data.py (normalize the estimated parameters whose cost function is greater than C_thres).
	2. fast_density_cluster.py (Decision graph to compute clusters according to Density Peak Algorithm 
	   in Ref. [2]) 




---------
References:
----------

Ref[1]: Thushara P. Abeyweera, Ernesto Merino, Morgan Huse. Inhibitory signaling blocks activating receptor clustering and induces cytoskeletal retraction in natural killer cells. J Cell Biol, 192 (4): 675???690, (2011).

Ref.[2]: Rodriguez A, Laio A. Clustering by fast search and find of density peaks. science. 2014;344(6191):1492-6


