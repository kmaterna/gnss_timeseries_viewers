# GNSS Time Series Viewers

## Overview:
This library contains a set of tools to read GNSS time series and velocities, remove earthquake/antenna offsets, solve for slopes, remove seasonal terms using several algorithms, and plot time series, stacks, and maps.  It can read data from the University of Nevada Reno (.tenv3) format, the Plate Boundary Observatory / Network Of The Americas (.pos) format, and the USGS format (.rneu).


## Dependencies and Installation
For full functionality, this library requires the Python packages listed in requirements.txt. 
You can create a new Conda environment using requirements.txt, or you can manually ensure that you have each of those libraries installed.
If you're manually creating the environment, don't forget the pip libraries: [Tectonic-Utils](https://github.com/kmaterna/Tectonic_Utils), 
and the [Earthscope CLI](https://gitlab.com/earthscope/public/earthscope-cli) for downloading Earthscope time series after 2023. 
For the mapping utility, I use the [pygmt](https://github.com/GenericMappingTools/pygmt) library.

Git clone this library onto your local machine and add the parent directory storing "GNSS_TimeSeries_Viewers" on your $PYTHONPATH. 


## Accessing Prerequisite GNSS Data
This library requires a local copy of GNSS time series, velocities, steps, offsets, and/or hydrological loading models from 
at least one GNSS database in their online-provided formats. Each one is slightly different.  The file ```getting_gnss_data/metadata.txt``` describes some of how to locate and download these files.
Remember that Earthscope now requires a log-in and the Earthscope-CLI for wget access.  
 

The script ```getting_gnss_data/update_data_holdings.sh``` now automates the download of many repositories locally, which might take 
a few hours (and ~10Gb) if you're downloading every time series from every processing center.  See specific scripts for individual databases.  

## Configuration

So, you've got some GNSS data onto your system. Now, to organize it and pass it into the library. An example file structure after running several scripts above might look like this:

```commandline 
├── Earthscope_Data
│   ├── Offsets/
│   ├── Earthscope_Event_Files/
│   ├── Time_Series/
│   ├── Velocities/
│   ├── cwu_config.txt
│   ├── nmt_config.txt
│   └── pbo_config.txt
├── UNR_Data
│   ├── Offsets/
│   ├── Time_Series/
│   ├── Velocities/
│   ├── UNR_coords_current.txt
│   └── unr_config.txt
├── USGS_Data
│   ├── Metadata/
│   ├── Offsets/
│   ├── Time_Series/
│   ├── Velocities/
│   └── usgs_config.txt
└── config.txt
```
With these directories on your system, you need to create a **VERY IMPORTANT** configuration file ```config.txt``` that tells the library where to find each GNSS data source.
In this case, ```config.txt``` contains something like the following:

```commandline 
# Example Config file for velocities, time series, and offsets
#############################################
# 
[py-config]

gps_data_dir = /Users/Path/to/../GPS_POS_Data/
pbo_config = /Users/Path/to/../GPS_POS_DATA/Earthscope_Data/pbo_config.txt
cwu_config = /Users/Path/to/../GPS_POS_DATA/Earthscope_Data/cwu_config.txt
nmt_config = /Users/Path/to/../GPS_POS_DATA/Earthscope_Data/nmt_config.txt
unr_config = /Users/Path/to/../GPS_POS_DATA/UNR_Data/unr_config.txt
usgs_config = /Users/Path/to/../GPS_POS_DATA/USGS_Data/usgs_config.txt
hydro_config = /Users/Path/to/../GPS_POS_DATA/Hydro/hydro_config.txt
```

and each of the inner config files reflect what is inside each database (see templates in the ```examples/``` directory). 
Please follow a similar templates for your own system. The path to the ```config.txt``` file will be passed into the library each time you use it.


## Library Features:
From this library, one can:
* Download local copies of GNSS velocities, time series, and offsets from the UNR, Earthscope, and USGS databases.
* Plot single GNSS time series (3-component)
* Read time series and velocities from UNR, PBO, CWU, and NMT formats
* Extract offsets at the times of earthquakes and antenna changes
* Remove seasonal terms using:
    * Least Squares fitting to a model with sines and cosines at 12-month and 6-month periods
    * A notch filter with notches at 12-month and 6-month periods
    * NLDAS, a high-resolution soil moisture and hydrological loading model for North America provided by UNAVCO
    * GLDAS, a medium-resolution soil moisture and hydrological loading model for North America provided by UNAVCO
    * LSDM, a high-resolution surface water and hydrological loading model globally provided by GFZ
* Plot stacks of time series
* Manually exclude bad stations from stacks if they are on the user-defined blacklist
* Calculate the formal velocity uncertainty and an empirical uncertainty on velocity estimates using the Allen Variance of the Rate (Hackl et al., 2011)
* Write out processed time series as text files in .pos format
* Create more complex experiments with GNSS time series (custom offsets, custom stacks, etc.)
* Specify the configuration of your system and the location of the local GNSS data with a plain-text config file, passed into the library each time it is used.

## Contributing
If you're using this library, I encourage you to reach out with feedback!  This includes bugs you find and features you'd like to see included. 
I'm happy to discuss and work together on this code. To contribute, you can submit issues or pull requests.


## Examples

### Example 1: How many stations within a region?
Using one function, we can figure out which stations in the PBO newtwork are within a certain radius or region. In this example, we read the PBO stations and return any within 100 km of a chosen coordinate.
```bash
example $ example_driver_find_within_radius.py 
Reading /Users/path/to/Velocity_Files/NAM08_pbovelfile_feb2018.txt
Returning 26 stations that are within 100.000 km of center -122.0000, 40.0000
['ORVB', 'P270', 'P272', 'P333', 'P334', 'P335', 'P336', 'P337', 'P339', 'P341', 'P344', 'P345', 'P346', 
'P349', 'P664', 'P665', 'P666', 'P667', 'P668', 'P669', 'P670', 'P671', 'P794', 'QUIN', 'SUTB', 'WDCB']
```


### Example 2: Single GNSS Time Series
As another example, we plot the PBO time series of P511, a station in Southern California with seasonal terms and antenna changes removed but earthquakes still left in. One possible usage of one driver is illustrated and its command-line outputs are shown.
```bash
example $ example_driver_single_plot.py P511
------- P511 --------
Viewing station P511, earthquakes_remove=0, outliers_remove=1, seasonals_remove=1, datasource=cwu, refframe=NA

Station P511: 
Reading /Users/path/to/Earthscope_Data/P511.cwu.final_nam14.pos
Reading data for station P511 in time range 2005-06-30:2019-11-30
Offset table for station P511:
 P511  2010 04 04  22 40        -7.66     0.52        0.55     0.40       3.19     1.87  OffEq ! EQ GU Location   32.14485  244.62646 ID ANSS(ComCat) ci14607652 


Earthquake table for station P511:
/Users/path/to/Earthscope_Event_Files/cwu_100404_2240_eqgu_coseis_kalts.evt: 244.70390  33.88694     0.55    -7.66      0.40     0.52   0.000     3.19     1.87  P511_GGU

Removing seasonals by lssq method.
Saving figure as P511_noseasons_lssq_cwu_NA_ts.jpg 
```
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/P511_noseasons_lssq_cwu_NA_ts.jpg)

### Example 3: A GNSS Stack
We can also plot a list of GNSS time series as a stack.  An example is shown below for stations in the northern San Francisco Bay Area. The parameters can be set in ```example_driver_stack.py```.  You can select which kinds of plots you want, whether it's vertical, horizontal, filtered, earthquakes removed, etc.
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/NBay_-122.0_38.0_40_TS_noeq.png)

If you have the pygmt library installed, you can make a simple map of the stations in your region:
![map](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/NBay_-122.0_38.0_40_map.png)

