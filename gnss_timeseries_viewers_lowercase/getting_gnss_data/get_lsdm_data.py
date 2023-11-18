#!/usr/bin/env python
"""
This Python script takes a list of stations
And downloads the loading products for each station from the German loading website
Extraction script is here:
http://rz-vm115.gfz-potsdam.de:8080/repository/entry/show?entryid=362f8705-4b87-48d1-9d86-2cfd1a2b6ac9
I use the Center of Figure "HYDL" products.
Download the executable at the link, and place it in the same directory as this script.
Put your desired stations into a text file with row format:
 name longitude latitude startdate[YYYY-MM-DD] enddate[YYYY-MM-DD]
 Example Rows:
     ASHL -122.670169 42.180690 2007-04-03 2022-10-22
     BLY1 -121.049065 42.406845 2011-11-15 2022-10-12
     CABL -124.563347 42.836097 1997-08-18 2022-10-22
Set the config parameters at the bottom of this script for your desired run.
Then run!
Makes a lot of command-line noise and takes 20 seconds per station, but it works!
Script by K. Materna, April 21, 2019
"""

import subprocess
import datetime as dt


def get_stations(input_file, product, output_dir, log_file):
    """ Function to download a bunch of LSDM loading products from the GFZ website! """
    ifile = open(input_file, 'r');
    ofile = open(log_file, 'w');
    ofile.close();  # just to start the log file
    for line in ifile:
        ofile = open(log_file, 'a');  # adding to the log each time
        temp = line.split();
        if len(line.split()) == 0:
            continue;
        station_name, lon, lat = temp[0], temp[1], temp[2];
        start_date, end_date = dt.datetime.strptime(temp[3], "%Y-%m-%d"), dt.datetime.strptime(temp[4], '%Y-%m-%d');
        file_destination = output_dir + station_name+"_LSDM_hydro";  # .txt extension will be added automatically

        print("Getting LSDM Hydro for station %s " % station_name);
        print("./extractlatlon_bilinintp_remote "+product+" CF "+dt.datetime.strftime(start_date, "%Y-%m-%d")+" " +
              dt.datetime.strftime(end_date, "%Y-%m-%d")+" "+lat+" "+lon+" -o " +
              file_destination);
        subprocess.call(["./extractlatlon_bilinintp_remote "+product+" CF "+dt.datetime.strftime(start_date, "%Y-%m-%d")
                         + " "+dt.datetime.strftime(end_date, "%Y-%m-%d")+" "+lat+" "+lon+" -o " +
                         file_destination], shell=True);
        ofile.write("\n"+line.split('\n')[0] + " " + file_destination + ".txt Downloaded at " + str(dt.datetime.now()));
        ofile.close();
    ifile.close();
    return;


if __name__ == "__main__":
    input_file = "CA_OR.txt";  # set your input file with station specs
    output_dir = "../LSDM/";  # set where you want to put your LSDM loads.  Should end with '/'.
    product = "HYDL";  # OPTIONS: HYDL, NTAL, NTOL+NTAL, etc.
    get_stations(input_file, product, output_dir, log_file='log.txt');
