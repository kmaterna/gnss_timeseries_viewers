# GNSS Time Series Viewers

## Overview: 
This library contains a set of tools to read GNSS time series and velocities, remove earthquake/antenna offsets, solve for slopes, remove seasonal terms using several algorithms, and plot time series, stacks, and maps.  It can read data from the University of Nevada Reno (.tenv3) format, the Plate Boundary Observatory / Network Of The Americas (.pos) format, and the USGS format (.rneu).  


## Dependencies and Installation
This library requires basic Python packages such as numpy and matplotlib. For full downloading functions, it requires a few other pip Python dependencies: pandas, html5lib, beautifulsoup4. For the "stations_within_radius" function, it also requires a file called 'haversine' that is located in my Utilities library (https://github.com/kmaterna/Tectonic_Utils), so please install that as welll (can be done through pip). For the mapping utility, I use the pygmt library based on GMT-6 (https://github.com/GenericMappingTools/pygmt). It's highly optional, but you can make an automatic map if you have it.  

You can clone this library locallly and put the parent directory storing "GNSS_TimeSeries_Viewers" on your $PYTHONPATH. The drivers in examples_and_configs provide some examples on how to access the major functions in this library.  I expect this library will work on Mac and Linux systems, but it hasn't been tested on Windows. 


## Accessing Prerequisite GNSS Data
This library requires a local copy of GNSS time series, velocities, steps, and/or hydrological loading models from various online repositories in their online-provided formats. The "getting_gnss_data/metadata.txt" file describes some of how to locate and download these files. The script 'getting_gnss_data/update_data_holdings.sh' now automates the download of all of the repositories locally, which might take a few hours (and ~10Gb) if you're downloading every time series from every processing center. 

Once you have a local copy of desired GNSS products, you must create an overall ```config.txt``` that tells the library where all the databases are located on your system.  
I create one overall config file (e.g., ```config.txt```, with absolute paths) for all databases, and one specific config file (e.g., ```unr_config.txt``` with relative paths) for each database that I've downloaded. 
Examples for ```config.txt``` and ```unr_config.txt``` etc. are provided in the "examples_and_configs"; please follow a similar templates for your own system. The path to the ```config.txt``` file will be passed into the library each time you use it.  


## Contributing
If you're using this library, I encourage you to reach out with feedback!  This includes bugs you find and features you'd like to see included. I'm happy to discuss and work together on this code. To contribute, you can submit issues or pull requests, or you can reach out by email; my email is kmaterna [at] usgs [dot] gov.


## Library Features:
From this library, one can:
* Download local copies of GNSS velocities, time series, and offsets from the UNR, GAGE, and USGS databases. 
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


## Examples

### Example 1: How many stations within a region? 
Using one function, we can figure out which stations in the PBO newtwork are within a certain radius or region. In this example, we read the PBO stations and return any within 100 km of a chosen coordinate. 
```bash
(pygmt) MacBook-Pro-2:example kmaterna$ example_driver_find_within_radius.py 
Reading /Users/kmaterna/Documents/GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt
Returning 26 stations that are within 100.000 km of center -122.0000, 40.0000
['ORVB', 'P270', 'P272', 'P333', 'P334', 'P335', 'P336', 'P337', 'P339', 'P341', 'P344', 'P345', 'P346', 
'P349', 'P664', 'P665', 'P666', 'P667', 'P668', 'P669', 'P670', 'P671', 'P794', 'QUIN', 'SUTB', 'WDCB']
```


### Example 2: Single GNSS Time Series
As another example, we plot the PBO time series of P511, a station in Southern California with seasonal terms and antenna changes removed but earthquakes still left in. One possible usage of one driver is illustrated and its command-line outputs are shown. 
```bash
(pygmt) MacBook-Pro-2:example kmaterna$ example_driver_single_plot.py P511
------- P511 --------
Viewing station P511, earthquakes_remove=0, outliers_remove=1, seasonals_remove=1, datasource=cwu, refframe=NA

Station P511: 
Reading /Users/kmaterna/Documents/GPS_POS_DATA/PBO_Data/P511.cwu.final_nam14.pos
Reading data for station P511 in time range 2005-06-30:2019-11-30
Offset table for station P511:
 P511  2010 04 04  22 40        -7.66     0.52        0.55     0.40       3.19     1.87  OffEq ! EQ GU Location   32.14485  244.62646 ID ANSS(ComCat) ci14607652 


Earthquake table for station P511:
/Users/kmaterna/Documents/GPS_POS_DATA/PBO_Event_Files/cwu_100404_2240_eqgu_coseis_kalts.evt: 244.70390  33.88694     0.55    -7.66      0.40     0.52   0.000     3.19     1.87  P511_GGU

Removing seasonals by lssq method.
Saving figure as P511_noseasons_lssq_cwu_NA_ts.jpg 
```
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/P511_noseasons_lssq_cwu_NA_ts.jpg)

### Example 3: A GNSS Stack
We can also plot a list of GNSS time series as a stack.  An example is shown below for stations in the northern San Francisco Bay Area. The parameters can be set in the 'driver_stack.py' file.  In the gps_stack.py driver, you can select which kinds of plots you want, whether it's vertical, horizontal, filtered, earthquakes removed, etc.
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/NBay_-122.0_38.0_40_TS_noeq.png)

If you have the pygmt library installed, you can make a simple map of the stations in your region: 
![map](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/examples_and_configs/example_pngs/NBay_-122.0_38.0_40_map.png)


