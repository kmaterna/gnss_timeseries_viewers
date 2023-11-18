#!/bin/bash
# Updating all data holdings as of June 2023.
# for USGS, be in an environment that has pandas
# for CWU, be in an environment that has "pip install earthscope-cli"
# Call this script from the parent directory where GNSS data is locally stored (ex: GPS_POS_DATA/)

# Setup
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# UPDATE EARTHSCOPE DATA FROM CWU PROCESSING CENTER
./get_earthscope_data.sh


# UPDATE UNR DATA: OFFSETS, COORDINATES, AND TIME SERIES
# VELOCITIES get updated manually with the datestamp and #URL placed on the first line.
unr_coords_cache_file="UNR_coords.txt"   # update this with the datestamp if desired
mkdir -p UNR_Data/
mkdir -p UNR_Data/Offsets/
cd UNR_Data/Offsets/
wget http://geodesy.unr.edu/NGLStationPages/steps.txt -O UNR_steps.txt
cd ../../

cd UNR_Data/
wget http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt -O $unr_coords_cache_file
cd ../

# Get time series from the stations in cache file.
mkdir -p UNR_Data/Time_Series/
cd UNR_Data/Time_Series/
python $DIR/get_unr_time_series.py ../$unr_coords_cache_file
cd ../../



# UPDATE USGS DATA: TIME SERIES AND VELOCITIES
mkdir -p USGS_Data/
cd USGS_Data/
python $DIR/get_usgs_data.py
cd ../



# UPDATE PBO HYDROLOGICAL MODELS
mkdir -p Hydro
mkdir -p Hydro/NLDAS/
cd Hydro/NLDAS/
wget -N --recursive --no-directories --no-parent --accept "*hyd" -e robots=off --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/hydro/nldas2/
cd ../../

mkdir -p Hydro/GLDAS/
cd Hydro/GLDAS/
wget -N --recursive --no-directories --no-parent --accept "*hyd" -e robots=off --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/hydro/gldas2/
cd ../../

mkdir -p Hydro/NOAH025
cd Hydro/NOAH025/
wget -N --recursive --no-directories --no-parent --accept "*hyd" -e robots=off --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/hydro/noah025/
cd ../../

