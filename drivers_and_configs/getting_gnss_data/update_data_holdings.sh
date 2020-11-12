#!/bin/bash
# Updating my data holdings as of October 6, 2020
# First, get into the major directory for GNSS Data.
cd ../../GPS_POS_DATA/

# UPDATE CWU/PBO POS DATA
cd PBO_Data/Time_Series/
wget -N --recursive --no-parent --no-directories --accept "*.cwu.final_igs14.pos, *.cwu.final_nam14.pos" ftp://data-out.unavco.org/pub/products/position
cd ../../

cd PBO_Data/PBO_Event_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu*coseis_kalts.evt, pbo*coseis_kalts.evt" ftp://data-out.unavco.org/pub/products/event/
cd ../../

cd PBO_Data/Offsets/
wget -N --recursive --no-parent --no-directories --accept "cwu*nam14.off, pbo*nam08.off" ftp://data-out.unavco.org/pub/products/offset/
cd ../../

cd PBO_Data/Velocity_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu.final_nam14.vel, cwu.final_igs14.vel" ftp://data-out.unavco.org/pub/products/velocity/
cd ../../



# UPDATE UNR DATA: OFFSETS, COORDINATES, AND TIME SERIES
cd UNR_Data/Offsets/
wget http://geodesy.unr.edu/NGLStationPages/steps.txt -O UNR_steps.txt
cd ../../

cd UNR_Data/
wget http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt -O UNR_coords_dec2020.txt
cd ../

# Check this script for parameter values and location on file system
cd UNR_Data/Time_Series/
python ../../../drivers_and_configs/getting_gnss_data/get_unr_time_series.py 
cd ../../



# UPDATE USGS DATA: TIME SERIES AND VELOCITIES
cd USGS_Data/
python ../../drivers_and_configs/getting_gnss_data/get_usgs_data.py
cd ../



# UPDATE PBO HYDROLOGICAL MODELS
cd Hydro/NLDAS/
wget -N --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/nldas2/
cd ../../

cd Hydro/GLDAS/
wget --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/gldas2
cd ../../

cd Hydro/NOAH025/
wget --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/noah025
cd ../../



