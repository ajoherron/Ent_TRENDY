# Instructions for running scripts
0. Set paths
1. 1_scaleacc_gij.sh
   1. Set the 'data_dir' path to where you've copied the raw data that needs processing.
   2. Set the 'gij_dir' path to an empty directory you've created to store all the scaled gij files.
   3. Run: 
      1. './1_scaleacc_gij.sh' 
2. 2_concatenate_time.py
   1. At the top of the script, update the 'GIJ_DIR' variable to the path you set for 'gij_dir' above.
   2. Similarly, update the 'INT_DIR' to a directory where you'd like to keep your intermediate NetCDFs.
   3. Run:
      1. python 2_concatenate_time.py
3. 3_concatenate_PFT.py
   1. 

# Description	


# Structure (from Nancy)
Subdirectory:
trendy-gcb2024/
  README
	Modules to load
	How-to run etc.
