#!/bin/bash
# Updating all data holdings as of June 2023.
# for USGS, be in an environment that has pandas
# for CWU, be in an environment that has "pip install earthscope-cli"
# Call this script from the parent directory where GNSS data is locally stored (ex: GPS_POS_DATA/)

# Setup
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


# UPDATE CWU/PBO POS, OFFSET, EVENT, AND VELOCITY DATA
# First, enable your es sso token with 'es sso login' (expires every day?)
# THIS WILL TAKE A WHILE, and it will look like it's not working at first, but it will be working.
# Breadth first search
mkdir -p PBO_Data
mkdir -p PBO_Data/Time_Series/
cd PBO_Data/Time_Series/
wget --recursive --no-directories --no-parent -N --accept="*.cwu.final_nam14.pos,*.cwu.final_igs14.pos" -e robots=off --level=2 --verbose --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/position/
cd ../../

mkdir -p PBO_Data/PBO_Event_Files/
cd PBO_Data/PBO_Event_Files/
wget -r -np --reject=tmp,ps,index* --accept "*_kalts.evt" -e robots=off --no-directories --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/event/
cp data.unavco.org/archive/gnss/products/event/*.evt .
rm -r data.unavco.org
cd ../../

mkdir -p PBO_Data/Offsets/
cd PBO_Data/Offsets/
wget -N --recursive --no-parent --no-directories --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/offset/cwu.kalts_nam14.off
cd ../../

mkdir -p PBO_Data/Velocities/
cd PBO_Data/Velocities/
wget -N --recursive --no-parent -e robots=off --no-directories --accept "cwu.final_nam14.vel, cwu.final_igs14.vel" --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/velocity/
cd ../../


# UPDATE UNR DATA: OFFSETS, COORDINATES, AND TIME SERIES
# VELOCITIES get updated manually with the datestamp and #URL placed on the first line.
unr_coords_cache_file="UNR_coords_jun2023.txt"   # update this with the datestamp
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

