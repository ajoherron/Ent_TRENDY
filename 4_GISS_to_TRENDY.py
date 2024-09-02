# Standard library imports
import xarray as xr
import pandas as pd
import numpy as np
from tqdm import tqdm

# Local imports
from config import INT_DIR, OUT_DIR, TOPO_PATH, SOIL_PATH, GIJ_PATH, PFT_PATH

#############################################
### Read in data as datasets using xarray ###
#############################################

# Standard variables
ds = xr.open_dataset(f"{INT_DIR}ds_standard_vars.nc")

# Topographic variables
ds_topo = xr.open_dataset(TOPO_PATH)

# Soil variables
ds_soil = xr.open_dataset(SOIL_PATH)

# Miscellaneous GIJ file used for axyp
ds_gij = xr.open_dataset(GIJ_PATH)

# File used for PFT naming conventions
ds_PFT_names = xr.open_dataset(PFT_PATH)

# Open Ent variable datasets
ds_ra001 = xr.open_dataset(f"{INT_DIR}ra001_PFT_dim.nc")
ds_ra003 = xr.open_dataset(f"{INT_DIR}ra003_PFT_dim.nc")
ds_ra005 = xr.open_dataset(f"{INT_DIR}ra005_PFT_dim.nc")
ds_ra007 = xr.open_dataset(f"{INT_DIR}ra007_PFT_dim.nc")
ds_ra016 = xr.open_dataset(f"{INT_DIR}ra016_PFT_dim.nc")
ds_ra017 = xr.open_dataset(f"{INT_DIR}ra017_PFT_dim.nc")
ds_ra019 = xr.open_dataset(f"{INT_DIR}ra019_PFT_dim.nc")
ds_ra020 = xr.open_dataset(f"{INT_DIR}ra020_PFT_dim.nc")
ds_ra021 = xr.open_dataset(f"{INT_DIR}ra021_PFT_dim.nc")
ds_ra022 = xr.open_dataset(f"{INT_DIR}ra022_PFT_dim.nc")
ds_ra023 = xr.open_dataset(f"{INT_DIR}ra023_PFT_dim.nc")
ds_ra024 = xr.open_dataset(f"{INT_DIR}ra024_PFT_dim.nc")
ds_ra025 = xr.open_dataset(f"{INT_DIR}ra025_PFT_dim.nc")
ds_ra028 = xr.open_dataset(f"{INT_DIR}ra028_PFT_dim.nc")
ds_ra031 = xr.open_dataset(f"{INT_DIR}ra031_PFT_dim.nc")
ds_ra032 = xr.open_dataset(f"{INT_DIR}ra032_PFT_dim.nc")
ds_ra033 = xr.open_dataset(f"{INT_DIR}ra033_PFT_dim.nc")
ds_ra034 = xr.open_dataset(f"{INT_DIR}ra034_PFT_dim.nc")
ds_ra035 = xr.open_dataset(f"{INT_DIR}ra035_PFT_dim.nc")
ds_ra036 = xr.open_dataset(f"{INT_DIR}ra036_PFT_dim.nc")
ds_ra037 = xr.open_dataset(f"{INT_DIR}ra037_PFT_dim.nc")
ds_ra038 = xr.open_dataset(f"{INT_DIR}ra038_PFT_dim.nc")
ds_ra039 = xr.open_dataset(f"{INT_DIR}ra039_PFT_dim.nc")
ds_ra040 = xr.open_dataset(f"{INT_DIR}ra040_PFT_dim.nc")

# Convert datasets to data arrays
ra001 = ds_ra001["ra001"]
ra003 = ds_ra003["ra003"]
ra005 = ds_ra005["ra005"]
ra007 = ds_ra007["ra007"]
ra016 = ds_ra016["ra016"]
ra017 = ds_ra017["ra017"]
ra019 = ds_ra019["ra019"]
ra020 = ds_ra020["ra020"]
ra021 = ds_ra021["ra021"]
ra022 = ds_ra022["ra022"]
ra023 = ds_ra023["ra023"]
ra024 = ds_ra024["ra024"]
ra025 = ds_ra025["ra025"]
ra028 = ds_ra028["ra028"]
ra031 = ds_ra031["ra031"]
ra032 = ds_ra032["ra032"]
ra033 = ds_ra033["ra033"]
ra034 = ds_ra034["ra034"]
ra035 = ds_ra035["ra035"]
ra036 = ds_ra036["ra036"]
ra037 = ds_ra037["ra037"]
ra038 = ds_ra038["ra038"]
ra039 = ds_ra039["ra039"]
ra040 = ds_ra040["ra040"]

#################
### Constants ###
#################

# ModelE does not account for leap years
days_in_year = 365

# Variable based on days in given month
# (Should work with actual data with varying months,
# all currently have January as default)
days_in_month = xr.DataArray(
    pd.to_datetime(ds["time"].values).days_in_month,
    coords=[ds["time"]],
    name="days_in_month",
)

######################################
### Formatting / Metadata / Saving ###
######################################


def replace_zeros_with_nan(ds):
    # If the input is a DataArray, handle it directly
    if isinstance(ds, xr.DataArray):
        return ds.where(ds != 0, np.nan)

    # If the input is a Dataset, create a copy and modify it
    ds_copy = ds.copy()

    # Loop through all data variables in the dataset
    for var in ds_copy.data_vars:
        # Replace 0s with NaNs
        ds_copy[var] = ds_copy[var].where(ds_copy[var] != 0, np.nan)

    return ds_copy


def format_coordinates_metadata(ds):
    ds_renamed = ds.copy()

    # Check if lon and lat are coordinates and rename if necessary
    if "lon" in ds_renamed.coords:
        ds_renamed = ds_renamed.rename({"lon": "longitude"})

    if "lat" in ds_renamed.coords:
        ds_renamed = ds_renamed.rename({"lat": "latitude"})

    # Check coordinates, reorder accordingly
    if len(ds_renamed.coords) == 3:
        ds_renamed = ds_renamed.transpose("longitude", "latitude", "time")
    elif len(ds_renamed.coords) == 4:
        if "Pool" in ds_renamed.coords:
            ds_renamed = ds_renamed.transpose("longitude", "latitude", "Pool", "time")
        elif "stlayer" in ds_renamed.coords:
            ds_renamed = ds_renamed.transpose(
                "longitude", "latitude", "stlayer", "time"
            )
        elif "smlayer" in ds_renamed.coords:
            ds_renamed = ds_renamed.transpose(
                "longitude", "latitude", "smlayer", "time"
            )
        elif "PFT" in ds_renamed.coords:
            ds_renamed = ds_renamed.transpose("longitude", "latitude", "PFT", "time")

    # Ensure metadata is global, rather than associated with individual variables
    ds_renamed = ds_renamed.to_dataset()

    # Update metadata
    ds_renamed.attrs["model_description"] = (
        "Ent Terrestrial Biosphere Model, configuration with NASA Goddard Institute for Space Studies ModelE-2.1 land surface hydrology physics (Kelley et al. 2020, doi:10.1029/2019MS002025), biophysics-only mode with vegetation structure and land use change as documented in Ito et al. (2020, doi:10.1029/2019MS002030), with the following updates to ModelE-2.1:  addition of a standalone land surface model (SLSM) driver for 0.5x0.5 degree simulations"
    )
    ds_renamed.attrs["institution"] = "NASA Goddard Institute for Space Studies"
    ds_renamed.attrs["contact"] = "Nancy.Y.Kiang@nasa.gov"

    return ds_renamed


def set_fill_values(ds):

    # Remove fill values for latitude/longitude
    ds.latitude.encoding["_FillValue"] = None
    ds.longitude.encoding["_FillValue"] = None

    # Set fill value for all variables to -99999
    for var in ds.data_vars:
        ds[var].encoding["_FillValue"] = -99999

    return ds


def rename_PFT(ds, ds_PFT_names):

    ### Rename PFT from numeric to categorical ###
    # Collect list of formatted netcdf names
    netcdf_name_list = []
    for netcdf_name in list(ds_PFT_names.lctype.values[:16]):
        cleaned_string = netcdf_name.decode("utf-8").strip()
        netcdf_name_list.append(cleaned_string)

    # Create netcdf name mapping dictionary, assign to PFT coordinates
    pft_netcdf_name_mapping = {}
    for i in range(len(netcdf_name_list)):
        pft_netcdf_name_mapping[str(i + 1)] = netcdf_name_list[i]
    ds = ds.assign_coords(
        PFT=[pft_netcdf_name_mapping[str(pft)] for pft in ds["PFT"].values]
    )

    ### Create PFT_longname data variable ###
    # Create a mapping between lctype and PFT
    lctype_to_pft = dict(zip(ds_PFT_names.lctype.values, ds.PFT.values))

    # Create a new DataArray with PFT as dimension instead of lctype
    new_longname = xr.DataArray(
        data=np.zeros(16, dtype="|S13"), dims=["PFT"], coords={"PFT": ds.PFT}
    )

    # Fill the new DataArray with lctype_longname values
    for lctype, pft in lctype_to_pft.items():
        new_longname.loc[pft] = ds_PFT_names.lctype_longname.sel(lctype=lctype).values

    # Decode the byte strings to regular strings
    decoded_longnames = [
        name.decode("utf-8").strip() if isinstance(name, bytes) else name.strip()
        for name in new_longname.values
    ]

    # Reassign the decoded longnames back to the DataArray using dtype=str to avoid converting back to byte strings
    ds["PFT_longname"] = xr.DataArray(
        data=np.array(
            decoded_longnames, dtype=str
        ),  # Use dtype=str to keep regular strings
        dims=["PFT"],
        coords={"PFT": ds.PFT},
    )
    return ds


def save_to_netcdf(ds_list):

    # Loop through all data arrays
    for ds in tqdm(ds_list):

        # Format coordinates into correct order and label
        ds = format_coordinates_metadata(ds)

        # Ensure all fill values are -99999
        ds = set_fill_values(ds)

        # Rename PFT coordinates and add longname variable
        if "PFT" in ds.coords:
            ds = rename_PFT(ds)

        # Save to output directory
        variable_name = list(ds.data_vars)[0]
        ds.to_netcdf(
            f"/discover/nobackup/aherron1/TRENDY/{OUT_DIR}/TEST_{variable_name}.nc",
            encoding={var: {"zlib": True, "complevel": 1} for var in ds.data_vars},
        )


########################
### ModelE Variables ###
########################

# Row 117
axyp = ds_gij["axyp"].copy()
axyp = axyp.rename("axyp")
axyp.attrs["units"] = "m^2"
axyp.attrs["long_name"] = "Grid cell area"

# Row 118
ocnfr = ds_topo["focean"].copy()
ocnfr = ocnfr.rename("ocnfr")
ocnfr.attrs["units"] = "fraction"
ocnfr.attrs["long_name"] = "Ocean fraction of grid cell"

# Row 119
soilfr = ds["soilfr"].copy()
soilfr /= 100
soilfr = soilfr.rename("soilfr")
soilfr.attrs["units"] = "fraction"
soilfr.attrs["long_name"] = "Soil including vegetation fraction of grid cell"

# Row 120
lakefr = ds_topo["flake"].copy()
lakefr = lakefr.rename("lakefr")
lakefr.attrs["units"] = "fraction"
lakefr.attrs["long_name"] = "Lake fraction of grid cell"

# Row 121
landicefr = ds_topo["fgice"].copy()
landicefr = landicefr.rename("landicefr")
landicefr.attrs["units"] = "fraction"
landicefr.attrs["long_name"] = "Permanent land ice fraction of grid cell"

# Row 122
bsfr = ds["bsfr"].copy()
bsfr /= 100
bsfr = bsfr.rename("bsfr")
bsfr.attrs["units"] = "fraction"
bsfr.attrs["long_name"] = "Bare component of land soil where bsfr+vsfr=soilfr"

# Row 123
vsfr = ds["vsfr"].copy()
vsfr /= 100
vsfr = vsfr.rename("vsfr")
vsfr.attrs["units"] = "fraction"
vsfr.attrs["long_name"] = "Vegetated component of land soil where bsfr+vsfr=soilfr"

# Row 124
dz = ds_soil["dz"].copy()
dz = dz.rename("dz")
dz.attrs["units"] = "m"
dz.attrs["long_name"] = "Thickness of each soil layer"

# Format and save each variable
variable_list = [
    axyp,
    ocnfr,
    soilfr,
    lakefr,
    landicefr,
    bsfr,
    vsfr,
    dz,
]
save_to_netcdf(variable_list)

##################
### Priority 1 ###
##################

# Row 12
tas = ds["tsurf"].copy()
tas += 273.15
tas = tas.rename("tas")
tas.attrs["units"] = "K"
tas.attrs["long_name"] = "Near-Surface Air Temperature"

# Row 13
pr = ds["prec"].copy()
pr /= 86400
pr = pr.rename("pr")
pr.attrs["units"] = "kg m-2 s-1"
pr.attrs["long_name"] = "Precipitation"

# Row 14
rsds = ds["incsw_grnd"].copy()
rsds = rsds.rename("rsds")
rsds.attrs["units"] = "W m-2"
rsds.attrs["long_name"] = "Surface Downwelling Shortwave Radiation"

# Row 15
mrso = ds["gwtr"].copy()
mrso = mrso.rename("mrso")
mrso.attrs["units"] = "kg m-2"
mrso.attrs["long_name"] = "Total Soil Moisture Content"

# Row 16
runoff_ugrnd = ds["runoff_ugrnd"].copy()
runoff_soil = ds["runoff_soil"].copy()
mrro = runoff_ugrnd + runoff_soil
mrro /= 86400
mrro = mrro.rename("mrro")
mrro.attrs["units"] = "kg m-2 s-1"
mrro.attrs["long_name"] = "Total Runoff"

# Row 17
evapotrans = ds["evap_land"].copy()
evapotrans = evapotrans * soilfr.copy() / (100 - ocnfr.copy()) / 86400
evapotrans = evapotrans.rename("evapotrans")
evapotrans.attrs["units"] = "kg m-2 s-1"
evapotrans.attrs["long_name"] = "Total Evapo-Transpiration"

# Row 21
transpft = ra005.copy() * 18.01528 * 0.001 * 2.5e6
transpft = transpft.rename("transpft")
transpft.attrs["units"] = "W m-2"
transpft.attrs["long_name"] = "Vegtype level transpiration"

# Format and save each variable
variable_list = [
    tas,
    pr,
    rsds,
    mrso,
    mrro,
    transpft,
    evapotrans,
    transpft,
]
save_to_netcdf(variable_list)

#####################################
### Land Carbon Variables (Pools) ###
#####################################

# Row 29
cVeg = ra001.copy() * (ra017.copy() - ra024.copy())
cVeg = cVeg.sum(dim="PFT")
cVeg = cVeg * vsfr.copy() / (100 - ocnfr.copy())
cVeg = cVeg.rename("cVeg")
cVeg.attrs["units"] = "kg m-2"
cVeg.attrs["long_name"] = "Carbon in Vegetation"

# Row 30
cLitter = ra001.copy() * (ra032.copy() + ra033.copy() + ra036.copy())
cLitter = cLitter.sum(dim="PFT")
cLitter = cLitter * vsfr.copy() / (100 - ocnfr.copy())
cLitter = cLitter.rename("cLitter")
cLitter.attrs["units"] = "kg m-2"
cLitter.attrs["long_name"] = "Carbon in Above-ground Litter Pool"

# Row 31
cSoil = ra001.copy() * (
    ra034.copy()
    + ra035.copy()
    + ra037.copy()
    + ra038.copy()
    + ra039.copy()
    + ra040.copy()
)
cSoil = cSoil.sum(dim="PFT")
cSoil = cSoil * vsfr.copy() / (100 - ocnfr.copy())
cSoil = cSoil.rename("cSoil")
cSoil.attrs["units"] = "kg m-2"
cSoil.attrs["long_name"] = "Carbon in Soil (including below-ground litter)"

# Row 32
ecvf = ds["ecvf"].copy()
cProduct = ecvf.copy() * 0.001 * days_in_year * 100 / (100 - ocnfr.copy())
cProduct = cProduct.rename("cProduct")
cProduct.attrs["units"] = "kg m-2"
cProduct.attrs["long_name"] = "Carbon in Products of Land Use Change"

# Row 33
cVegpft = ra017.copy() - ra024.copy()
cVegpft = cVegpft.rename("cVegpft")
cVegpft.attrs["units"] = "kg m-2"
cVegpft.attrs["long_name"] = "Vegtype level Carbon in Vegetation"

# Row 34
cSoilpft = (
    ra034.copy()
    + ra035.copy()
    + ra037.copy()
    + ra038.copy()
    + ra039.copy()
    + ra040.copy()
)
cSoilpft = cSoilpft.rename("cSoilpft")
cSoilpft.attrs["units"] = "kg m-2"
cSoilpft.attrs["long_name"] = "Vegtype level Carbon in Vegetation"

# Format and save each variable
variable_list = [
    cVeg,
    cLitter,
    cSoil,
    cProduct,
    cVegpft,
    cSoilpft,
]
save_to_netcdf(variable_list)

######################################
### Land Carbon Variables (Fluxes) ###
######################################

# Row 36
gpp = ds["gpp"].copy()
gpp = gpp * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
gpp = gpp.rename("gpp")
gpp.attrs["units"] = "kg m-2 s-1"
gpp.attrs["long_name"] = "Gross Primary Production"

# Row 37
rauto = ds["rauto"].copy()
ra = rauto * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
ra = ra.rename("ra")
ra.attrs["units"] = "kg m-2 s-1"
ra.attrs["long_name"] = "Autotrophic (Plant) Respiration"

# Row 38
npp = gpp.copy() - ra.copy()
npp = npp * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
npp = npp.rename("npp")
npp.attrs["units"] = "kg m-2 s-1"
npp.attrs["long_name"] = "Net Primary Production"

# Row 39
soilresp = ds["soilresp"].copy()
rh = soilresp * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
rh = rh.rename("rh")
rh.attrs["units"] = "kg m-2 s-1"
rh.attrs["long_name"] = "Heterotrophic Respiration"

# Row 42
ecvf = ds["ecvf"].copy()
fLuc = ecvf * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
fLuc = fLuc.rename("fLuc")
fLuc.attrs["units"] = "kg m-2 s-1"
fLuc.attrs["long_name"] = "CO2 Flux to Atmosphere from Land Use Change"

# Row 43
gpp = ds["gpp"].copy()
nbp = gpp.copy() - rauto.copy() - soilresp.copy() - ecvf.copy()
nbp = nbp * 0.001 / 86400 * soilfr.copy() / (100 - ocnfr.copy())
nbp = nbp.rename("nbp")
nbp.attrs["units"] = "kg m-2 s-1"
nbp.attrs["long_name"] = "Net Biospheric Production"

# Row 45
gpppft = ra003.copy() * 0.001 / 86400
gpppft = gpppft.rename("gpppft")
gpppft.attrs["units"] = "kg m-2 s-1"
gpppft.attrs["long_name"] = "Vegtype level GPP"

# Row 46
npppft = ra003.copy() - ra016.copy()
npppft = npppft * 0.001 / 86400
npppft = npppft.rename("npppft")
npppft.attrs["units"] = "kg m-2 s-1"
npppft.attrs["long_name"] = "Vegtype level NPP"

# Row 47
rhpft = ra025.copy() * 0.001 / 86400
rhpft = rhpft.rename("rhpft")
rhpft.attrs["units"] = "kg m-2 s-1"
rhpft.attrs["long_name"] = "Vegtype level Rh"

# Row 50
landCoverFrac = ra001.copy() * soilfr.copy() / (100 - ocnfr.copy())
landCoverFrac = landCoverFrac.rename("landCoverFrac")
landCoverFrac.attrs["units"] = "fraction"
landCoverFrac.attrs["long_name"] = "Fractional Land Cover of PFT"

# Row 51
oceanCoverFrac = ocnfr.copy()
oceanCoverFrac = oceanCoverFrac.rename("oceanCoverFrac")
oceanCoverFrac.attrs["units"] = "fraction"
oceanCoverFrac.attrs["long_name"] = "Fractional Ocean Cover"

# Row 53
lai = ds["LAI"].copy()
lai = lai.rename("lai")
lai.attrs["units"] = ""
lai.attrs["long_name"] = "Leaf Area Index"

# Row 54
laipft = ra007.copy()
laipft = laipft.rename("laipft")
laipft.attrs["units"] = ""
laipft.attrs["long_name"] = "Vegtype level Leaf Area Index"

# Row 57
cLeaf = ra001.copy() * ra019.copy()
cLeaf = cLeaf.sum(dim="PFT")
cLeaf = cLeaf * soilfr.copy() / (100 - ocnfr.copy())
cLeaf = cLeaf.rename("cLeaf")
cLeaf.attrs["units"] = "kg m-2"
cLeaf.attrs["long_name"] = "Carbon in Leaves"

# Row 58
cWood = ra020.copy() * ra021.copy()
cWood = cWood.sum(dim="PFT")
cWood = cWood * soilfr.copy() / (100 - ocnfr.copy())
cWood = cWood.rename("cWood")
cWood.attrs["units"] = "kg m-2"
cWood.attrs["long_name"] = "Carbon in Woood"

# Row 59
cRoot = ra022.copy() * ra023.copy()
cRoot = cRoot.sum(dim="PFT")
cRoot = cRoot * soilfr.copy() / (100 - ocnfr.copy())
cRoot = cRoot.rename("cRoot")
cRoot.attrs["units"] = "kg m-2"
cRoot.attrs["long_name"] = "Carbon in Roots"

# Row 60
cCwd = ra036.copy()
cCwd = cCwd.sum(dim="PFT")
cCwd = cCwd * soilfr.copy() / (100 - ocnfr.copy())
cCwd = cCwd.rename("cCwd")
cCwd.attrs["units"] = "kg m-2"
cCwd.attrs["long_name"] = "Carbon in Coarse Woody Debris"

# Row 62
fVegLitter = ra001 * ra031
fVegLitter = fVegLitter.sum(dim="PFT")
fVegLitter = fVegLitter * soilfr.copy() / (100 - ocnfr.copy())
fVegLitter = fVegLitter.rename("fVegLitter")
fVegLitter.attrs["units"] = "kg m-2 s-1"
fVegLitter.attrs["long_name"] = "Total Carbon Flux from Vegetation to Litter"

# Row 61

# Formula
C_soil_SURFMET = ra001.copy() * ra032.copy()
C_soil_SURFSTR = ra001.copy() * ra033.copy()
C_soil_SOILMET = ra001.copy() * ra034.copy()
C_soil_SOILSTR = ra001.copy() * ra035.copy()
C_soil_CWD = ra001.copy() * ra036.copy()
C_soil_SURFMIC = ra001.copy() * ra037.copy()
C_soil_SOILMIC = ra001.copy() * ra038.copy()
C_soil_SLOW = ra001.copy() * ra039.copy()
C_soil_PASSIVE = ra001.copy() * ra040.copy()

# Sum along PFT dimension
C_soil_SURFMET = C_soil_SURFMET.sum(dim="PFT")
C_soil_SURFSTR = C_soil_SURFSTR.sum(dim="PFT")
C_soil_SOILMET = C_soil_SOILMET.sum(dim="PFT")
C_soil_SOILSTR = C_soil_SOILSTR.sum(dim="PFT")
C_soil_CWD = C_soil_CWD.sum(dim="PFT")
C_soil_SURFMIC = C_soil_SURFMIC.sum(dim="PFT")
C_soil_SOILMIC = C_soil_SOILMIC.sum(dim="PFT")
C_soil_SLOW = C_soil_SLOW.sum(dim="PFT")
C_soil_PASSIVE = C_soil_PASSIVE.sum(dim="PFT")

# Unit convert multiplier
C_soil_SURFMET = C_soil_SURFMET * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SURFSTR = C_soil_SURFSTR * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SOILMET = C_soil_SOILMET * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SOILSTR = C_soil_SOILSTR * soilfr.copy() / (100 - ocnfr.copy())
C_soil_CWD = C_soil_CWD * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SURFMIC = C_soil_SURFMIC * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SOILMIC = C_soil_SOILMIC * soilfr.copy() / (100 - ocnfr.copy())
C_soil_SLOW = C_soil_SLOW * soilfr.copy() / (100 - ocnfr.copy())
C_soil_PASSIVE = C_soil_PASSIVE * soilfr.copy() / (100 - ocnfr.copy())

# Combine along Pool dimension
cSoilpools = xr.concat(
    [
        C_soil_SURFMET.rename("C_soil_SURFMET"),
        C_soil_SURFSTR.rename("C_soil_SURFSTR"),
        C_soil_SOILMET.rename("C_soil_SOILMET"),
        C_soil_SOILSTR.rename("C_soil_SOILSTR"),
        C_soil_CWD.rename("C_soil_CWD"),
        C_soil_SURFMIC.rename("C_soil_SURFMIC"),
        C_soil_SOILMIC.rename("C_soil_SOILMIC"),
        C_soil_SLOW.rename("C_soil_SLOW"),
        C_soil_PASSIVE.rename("C_soil_PASSIVE"),
    ],
    dim="Pool",
)
cSoilpools = cSoilpools.rename("cSoilpools")
cSoilpools["Pool"] = [
    "SURFMET",
    "SURFSTR",
    "SOILMET",
    "SOILSTR",
    "CWD",
    "SURFMIC",
    "SOILMIC",
    "SLOW",
    "PASSIVE",
]
cSoilpools.attrs["units"] = "kg m-2"
cSoilpools.attrs["long_name"] = "Carbon in Soil Pools"
cSoilpools = replace_zeros_with_nan(cSoilpools)

# Format and save each variable
variable_list = [
    gpp,
    gpppft,
    npppft,
    oceanCoverFrac,
    lai,
    ra,
    npp,
    rh,
    fLuc,
    nbp,
    rhpft,
    landCoverFrac,
    laipft,
    cLeaf,
    cWood,
    cRoot,
    cCwd,
    cSoilpools,
    fVegLitter,
]
save_to_netcdf(variable_list)

#######################
### Second Priority ###
#######################

# Row 78

# bs_tlay variables
bs_tlay1 = ds["bs_tlay1"].copy()
bs_tlay2 = ds["bs_tlay2"].copy()
bs_tlay3 = ds["bs_tlay3"].copy()
bs_tlay4 = ds["bs_tlay4"].copy()
bs_tlay5 = ds["bs_tlay5"].copy()
bs_tlay6 = ds["bs_tlay6"].copy()

# vs_tlay variables
vs_tlay1 = ds["vs_tlay1"].copy()
vs_tlay2 = ds["vs_tlay2"].copy()
vs_tlay3 = ds["vs_tlay3"].copy()
vs_tlay4 = ds["vs_tlay4"].copy()
vs_tlay5 = ds["vs_tlay5"].copy()
vs_tlay6 = ds["vs_tlay6"].copy()

# Formula & Unit convert multiplier
tsl_1 = (bsfr.copy() * bs_tlay1 + vsfr.copy() * vs_tlay1) / soilfr.copy() + 273.15
tsl_2 = (bsfr.copy() * bs_tlay2 + vsfr.copy() * vs_tlay2) / soilfr.copy() + 273.15
tsl_3 = (bsfr.copy() * bs_tlay3 + vsfr.copy() * vs_tlay3) / soilfr.copy() + 273.15
tsl_4 = (bsfr.copy() * bs_tlay4 + vsfr.copy() * vs_tlay4) / soilfr.copy() + 273.15
tsl_5 = (bsfr.copy() * bs_tlay5 + vsfr.copy() * vs_tlay5) / soilfr.copy() + 273.15
tsl_6 = (bsfr.copy() * bs_tlay6 + vsfr.copy() * vs_tlay6) / soilfr.copy() + 273.15

# Format
tsl = xr.concat([tsl_1, tsl_2, tsl_3, tsl_4, tsl_5, tsl_6], dim="stlayer")
tsl["stlayer"] = [1, 2, 3, 4, 5, 6]
tsl = tsl.rename("tsl")
tsl.attrs["units"] = "K"
tsl.attrs["long_name"] = "Temperature of Soil"

# Row 79

# bs_wlay variables
bs_wlay1 = ds["bs_wlay1"].copy()
bs_wlay2 = ds["bs_wlay2"].copy()
bs_wlay3 = ds["bs_wlay3"].copy()
bs_wlay4 = ds["bs_wlay4"].copy()
bs_wlay5 = ds["bs_wlay5"].copy()
bs_wlay6 = ds["bs_wlay6"].copy()

# Formula & Unit convert multiplier
msl_1 = (bsfr.copy() * bs_wlay1 + vsfr.copy() * vs_tlay1) / soilfr.copy()
msl_2 = (bsfr.copy() * bs_wlay2 + vsfr.copy() * vs_tlay2) / soilfr.copy()
msl_3 = (bsfr.copy() * bs_wlay3 + vsfr.copy() * vs_tlay3) / soilfr.copy()
msl_4 = (bsfr.copy() * bs_wlay4 + vsfr.copy() * vs_tlay4) / soilfr.copy()
msl_5 = (bsfr.copy() * bs_wlay5 + vsfr.copy() * vs_tlay5) / soilfr.copy()
msl_6 = (bsfr.copy() * bs_wlay6 + vsfr.copy() * vs_tlay6) / soilfr.copy()

# Format
msl = xr.concat([msl_1, msl_2, msl_3, msl_4, msl_5, msl_6], dim="smlayer")
msl["smlayer"] = [1, 2, 3, 4, 5, 6]
msl = msl.rename("msl")
msl.attrs["units"] = "kg m-2"
msl.attrs["long_name"] = "Moisture of Soil"

# Row 80
wetcan_evap = ds["wetcan_evap"].copy()
evspsblveg = wetcan_evap * vsfr.copy() / (100 - ocnfr.copy()) / 86400
evspsblveg = evspsblveg.rename("evspsblveg")
evspsblveg.attrs["units"] = "kg m-2 s-1"
evspsblveg.attrs["long_name"] = "Evaporation from Canopy"

# Row 81
evap_land = ds["evap_land"].copy()
drycan_evap = ds["drycan_evap"].copy()
evspsblsoi = (evap_land - (drycan_evap + wetcan_evap)) * vsfr.copy() / soilfr.copy()
evspsblsoi = evspsblsoi.copy() * soilfr.copy() / (100 - ocnfr.copy()) / 86400
evspsblsoi = evspsblsoi.rename("evspsblsoi")
evspsblsoi.attrs["units"] = "kg m-2 s-1"
evspsblsoi.attrs["long_name"] = "Water Evaporation from Soil"

# Row 82
drycan_evap = ds["drycan_evap"].copy()
tran = drycan_evap * vsfr.copy() / (100 - ocnfr.copy()) / 86400
tran = tran.rename("tran")
tran.attrs["units"] = "kg m-2 s-1"
tran.attrs["long_name"] = "Transpiration"

# Row 83
theightpft = ra003.copy().rename("theightpft")
theightpft.attrs["units"] = "m"
theightpft.attrs["long_name"] = "Vegtype level tree heights"

# Format and save each variable
variable_list = [
    tsl,
    msl,
    evspsblveg,
    evspsblsoi,
    tran,
    theightpft,
]
save_to_netcdf(variable_list)
