#!/bin/bash

# Set location of data that needs scaling
#data_dir="/discover/nobackup/aherron1/TRENDY/spinup2b"
data_dir="/discover/nobackup/aherron1/TRENDY/control_tr"

# Set where the scaled data should go
#gij_dir="/discover/nobackup/aherron1/TRENDY/spinup2b_gij"
gij_dir="/discover/nobackup/aherron1/TRENDY/control_tr_gij"

# Change to target directory
cd "$gij_dir" || exit

# Loop through all files in the data directory
for FILE in "$data_dir"/*; do
    if [ -f "$FILE" ]; then
        echo "Processing $FILE"
        scaleacc "$FILE" gij
    fi
done
