#!/bin/bash
# Updating my data holdings as of October 6, 2020
# First, get into the major directory for GNSS Data.
cd ../../GPS_POS_DATA/

# UPDATE CWU POS DATA
cd PBO_Data/
wget -N --recursive --no-parent --no-directories --accept "*.cwu.final_igs14.pos, *.cwu.final_nam14.pos" ftp://data-out.unavco.org/pub/products/position
cd ../

cd PBO_Event_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu*coseis_kalts.evt, pbo*coseis_kalts.evt" ftp://data-out.unavco.org/pub/products/event/
cd ../

cd Offsets/
wget -N --recursive --no-parent --no-directories --accept "cwu*nam14.off, pbo*nam08.off" ftp://data-out.unavco.org/pub/products/offset/
cd ../

cd Velocity_Files/
wget -N --recursive --no-parent --no-directories --accept "cwu.snips_nam14.vel, cwu.snaps_nam14.vel" ftp://data-out.unavco.org/pub/products/velocity/
cd ../

cd PBO_Hydro/NLDAS/
wget -N --recursive --no-parent --no-directories --accept "*.hyd" ftp://data-out.unavco.org/pub/products/hydro/nldas2/
cd ../../


# UPDATE UNR DATA: OFFSETS AND COORDINATES
cd Offsets/
wget http://geodesy.unr.edu/NGLStationPages/steps.txt -O UNR_steps.txt
cd ../

cd UNR_Data/
wget http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt -O UNR_coords_oct2020.txt
cd ../

# Check this script for parameter values
cd UNR_Data/
python get_unr_data.py 
cd ../
