#!/bin/bash
# Updating all data holdings as of October 2020.
# Call this script from the parent directory where GNSS data is locally stored (ex: GPS_POS_DATA/)

# Setup
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


# UPDATE CWU/PBO POS DATA
mkdir -p PBO_Data
mkdir -p PBO_Data/Time_Series/
cd PBO_Data/Time_Series/
wget -N --recursive --no-parent --no-directories --accept "*.cwu.final_igs14.pos, *.cwu.final_nam14.pos" ftp://data-out.unavco.org/pub/products/position
cd ../../

mkdir -p PBO_Data/PBO_Event_Files/
cd PBO_Data/PBO_Event_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu*coseis_kalts.evt, pbo*coseis_kalts.evt" ftp://data-out.unavco.org/pub/products/event/
cd ../../

mkdir -p PBO_Data/Offsets/
cd PBO_Data/Offsets/
wget -N --recursive --no-parent --no-directories --accept "cwu*nam14.off, pbo*nam08.off" ftp://data-out.unavco.org/pub/products/offset/
cd ../../

mkdir -p PBO_Data/Velocity_Files/
cd PBO_Data/Velocity_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu.final_nam14.vel, cwu.final_igs14.vel" ftp://data-out.unavco.org/pub/products/velocity/
cd ../../



# UPDATE UNR DATA: OFFSETS, COORDINATES, AND TIME SERIES
unr_coords_cache_file="UNR_coords_dec2020.txt"   # update this with the datestamp
mkdir -p UNR_Data/
mkdir -p UNR_Data/Offsets/
cd UNR_Data/Offsets/
wget http://geodesy.unr.edu/NGLStationPages/steps.txt -O UNR_steps.txt
cd ../../

cd UNR_Data/
wget http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt -O $unr_coords_cache_file
cd ../

# Get time series from the stations in cache file
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
wget -N --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/nldas2/
cd ../../

mkdir -p Hydro/GLDAS/
cd Hydro/GLDAS/
wget --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/gldas2
cd ../../

mkdir -p Hydro/NOAH025
cd Hydro/NOAH025/
wget --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/noah025
cd ../../



