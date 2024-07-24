######################
### User-set paths ###
######################

GIJ_DIR = "/discover/nobackup/aherron1/TRENDY/control_tr_gij/"
INT_DIR = "/discover/nobackup/aherron1/TRENDY/control_intermediates/"
OUT_DIR = "/discover/nobackup/aherron1/TRENDY/control_output/"

# Used ocnfr, lakefr, landicefr
TOPO_PATH = f"{INT_DIR}Z2HX2fromZ1QX1N.BS1.nc"

# Only used for dz
SOIL_PATH = f"{INT_DIR}S144X900098M.ext.nc"

# Can choose any GIJ file, as this is only used for axyp
GIJ_PATH = f"{INT_DIR}APR1701.gijL0_control_tr.nc"


#########################################
### No need to alter below this point ###
#########################################

STANDARD_VARIABLES = [
    "LAI",
    "tsurf",
    "prec",
    "incsw_grnd",
    "gwtr",
    "runoff_soil",
    "runoff_ugrnd",
    "evap_land",
    "ecvf",
    "gpp",
    "rauto",
    "soilresp",
    "soilfr",
    "wetcan_evap",
    "drycan_evap",
    "evap_land",
    "bsfr",
    "bs_tlay1",
    "bs_tlay2",
    "bs_tlay3",
    "bs_tlay4",
    "bs_tlay5",
    "bs_tlay6",
    "bs_wlay1",
    "bs_wlay2",
    "bs_wlay3",
    "bs_wlay4",
    "bs_wlay5",
    "bs_wlay6",
    "vsfr",
    "vs_tlay1",
    "vs_tlay2",
    "vs_tlay3",
    "vs_tlay4",
    "vs_tlay5",
    "vs_tlay6",
]


R_VARIABLES = [
    "ra001",
    "ra003",
    "ra005",
    "ra007",
    "ra016",
    "ra017",
    "ra019",
    "ra020",
    "ra021",
    "ra022",
    "ra023",
    "ra024",
    "ra025",
    "ra028",
    "ra031",
    "ra032",
    "ra033",
    "ra034",
    "ra035",
    "ra036",
    "ra037",
    "ra038",
    "ra039",
    "ra040",
]
