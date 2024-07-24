import xarray as xr
from tqdm import tqdm

from variables import R_VARIABLES
from paths import INT_DIR


# Concatenate ra___ variables along PFT dimension
def concat_PFT_variables(ENT, ds):

    variable_list = []
    variable_dict = {}

    for PFT in range(1, 17):

        # Create list of all variables
        if len(str(PFT)) == 1:
            variable = ENT + "00" + str(PFT)
        else:
            variable = ENT + "0" + str(PFT)
        variable_list.append(variable)

        # Create dictionary for renaming variables
        variable_dict[variable] = str(PFT)

    # Rename variables with dictionary
    ds_renamed = ds.rename(variable_dict).copy()

    # Create list of variable names
    new_var_names = list(variable_dict.values())

    # Concatenate the variables along a new dimension 'PFT'
    ds_pft = xr.concat([ds_renamed[var] for var in new_var_names], dim="PFT")

    # Add PFT coordinate
    ds_pft = ds_pft.assign_coords(PFT=("PFT", new_var_names))
    return ds_pft


# Loop through each variable, concatenate along PFT dimension for each, save as new dataset
for ra in tqdm(R_VARIABLES):
    ds_ra = xr.open_dataset(f"{INT_DIR}{ra}.nc")
    ra_concat_array = concat_PFT_variables(ra, ds_ra)
    ra_concat_dataset = ra_concat_array.to_dataset(name=ra)
    ra_concat_dataset.to_netcdf(f"{INT_DIR}{ra}_PFT_dim.nc")
