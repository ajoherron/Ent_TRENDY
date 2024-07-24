# Standard imports
import glob
import xarray as xr
import pandas as pd
from tqdm import tqdm

# Local imports
from config import GIJ_DIR, INT_DIR, STANDARD_VARIABLES, R_VARIABLES


months_to_mm = {
    "JAN": "01",
    "FEB": "02",
    "MAR": "03",
    "APR": "04",
    "MAY": "05",
    "JUN": "06",
    "JUL": "07",
    "AUG": "08",
    "SEP": "09",
    "OCT": "10",
    "NOV": "11",
    "DEC": "12",
}


# Expand dimensions along time
def add_time_dimension(ds, time_value):
    return ds.expand_dims(time=[time_value])


# Standardize fill/missing values to number specified by TRENDY protocol
def normalize_fill_values(ds, variables):
    for var in variables:
        if var in ds:
            ds[var].encoding["_FillValue"] = -99999
            ds[var].encoding["missing_value"] = -99999
    return ds


def concatenate_monthly_data(GIJ_DIR, INT_DIR, VARIABLES, file_name):

    # Gather all file paths
    GIJ_DIR += "*"
    file_paths = sorted(glob.glob(GIJ_DIR))
    print(len(file_paths))

    # Calculate adjusted time values
    time_values = []
    for file in tqdm(file_paths):

        # Parse month and year from file name
        month = file.split("/")[-1][:3]
        MM = months_to_mm[month]
        YYYY = file.split("/")[-1][3:7]

        # Add datetime value to list
        date_str = f"{YYYY}-{MM}-01"
        year_dt = pd.to_datetime(date_str, format="%Y-%m-%d", errors="coerce")
        time_values.append(year_dt)

    # Convert time_values to datetime64[D]
    time_values = pd.to_datetime(time_values).floor("D").values.astype("datetime64[D]")

    # Create subset dataset for each dataset in file path
    subset_datasets = []
    for file_path in tqdm(
        file_paths, desc="Selecting variables, normalizing fill values"
    ):
        ds = xr.open_dataset(file_path, decode_times=False)
        ds_subset = ds[VARIABLES]
        ds_subset = normalize_fill_values(ds_subset, VARIABLES)
        subset_datasets.append(ds_subset)

    # Add time dimension to all subset_datasets
    subset_datasets = [
        add_time_dimension(ds, time)
        for ds, time in tqdm(
            zip(subset_datasets, time_values), desc="Adding time dimension"
        )
    ]

    # Concatenate all subset datasets along time dimensions
    ds_all = xr.concat(
        subset_datasets,
        dim="time",
        coords="minimal",
        data_vars="minimal",
        compat="override",
    )

    # Sort dataset by time
    ds_all = ds_all.sortby("time")

    # Saved concatenated dataset in data folder
    ds_path = INT_DIR + f"{file_name}.nc"
    ds_all.to_netcdf(ds_path)
    print(f"Concatenated {len(file_paths)} datasets into:\n{ds_all}")


# Concatentate standard variables into a single dataset
file_name = "ds_standard_vars"
concatenate_monthly_data(GIJ_DIR, INT_DIR, STANDARD_VARIABLES, file_name)


# Create datasets for each Ent variable
# (including all PFT numbers for each)
for X in R_VARIABLES:

    # Loop through each PFT, add to variable list
    variable_subset = []
    for PFT in range(1, 17):
        if len(str(PFT)) == 1:
            variable = X + "00" + str(PFT)
        else:
            variable = X + "0" + str(PFT)
        variable_subset.append(variable)

    # Concatenate each Ent variable as a single dataset
    output_file_name = X
    concatenate_monthly_data(GIJ_DIR, INT_DIR, variable_subset, output_file_name)
