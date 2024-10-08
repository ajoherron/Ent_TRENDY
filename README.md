# Instructions for running scripts
0. Setup
   1. Set paths
      1. First, set all the relevant paths in the config file
      2. Next, set the paths for the input data directory and gij directory in 1_scaleacc_gij.sh
   2. Activate conda environment with 'conda env -f environment.yml'
      1. If you have access to the following dependencies, you can skip this step:
         1. Python
         2. glob2
         3. xarray
         4. pandas
         5. tqdm
1. ```1_scaleacc_gij.sh```
   1. This will loop though all the input files, and scale them using the following command:
      ```sh
      scaleacc <file> gij
2. ```2_concatenate_time.py```
   1. This will create NetCDF files for the relevant data:
      1. One file for standard variables (tas, prec, soilfr, etc.)
      2. One file for each Ent variable (ex: ra001)
   2. All these files will be concatenated along the time dimension so that they fit in a single file.
3. ```3_concatenate_PFT.py```
   1. This will concatenate each of the Ent variables along the PFT dimension. This will result in a NetCDF file with one variable and an added PFT coordinate, rather than having 16 separate variables.
4. ```4_GISS_to_TRENDY.py or 4_GISS_to_TRENDY.ipynb```
   1. This step has both python script and jupyter notebook options.
      1. The jupyter notebook option can be handy on JupyterHub, particularly when requesting a full node session for this step.
   2. This step converts the intermediate GISS variables created thus far into TRENDY variables, using formulae supplied. 
   3. As this is a data-intensive step, I've broken up the variables into chunk (same as those from the trendy_listofvariables_GCP2023 document).
   4. This will output the finished results in the directory specified in Step 0.
5. (Secondary) Remember to add README to output directory
   1. This was copied from the GlobalCarbonBudget-Protocol-2023 document.
   2. This has been labeled "TRENDY_README.md" to avoid confusion with the primary README for this repository.

# Structure (from Nancy)
Subdirectory:
-trendy-gcb2024/
	- README
		- Modules to load
		- How-to run etc.
