# GNSS Time Series Viewers

## Code Description: 
This library in GPS_TOOLS/ contains a set of tools to read GNSS time series and velocities, remove earthquake/antenna offsets, solve for slopes, remove seasonal terms using several algorithms, and plot time series, stacks, and maps.  It can read data from both the University of Nevada Reno (.tenv3) and the Plate Boundary Observatory / Network Of The Americas (.pos) formats. The tools here are meant to be modular, useful for stringing together into more complex experiments (custom offsets, custom stacks, etc.).  


## Features
* Plot single time teries (3-component)
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
* Write out processed time series as text files in .pos format


## Examples

### How many stations within a region? 
Using one function, we can figure out which stations are within a certain radius or region. In this example, we read the PBO stations and return any within 100 km of a chosen coordinate. 
```bash
(base) Kathryns-MacBook-Pro-2:example kmaterna$ driver_find_within_radius.py 
Reading /Users/kmaterna/Documents/GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt
Returning 26 stations that are within 100.000 km of center -122.0000, 40.0000
['ORVB', 'P270', 'P272', 'P333', 'P334', 'P335', 'P336', 'P337', 'P339', 'P341', 'P344', 'P345', 'P346', 'P349', 'P664', 'P665', 'P666', 'P667', 'P668', 'P669', 'P670', 'P671', 'P794', 'QUIN', 'SUTB', 'WDCB']
```


### Single GNSS Time Series
As one example, we plot the time series of P511, a station in Southern California with seasonal terms and antenna changes removed but earthquakes still left in. One possible usage of one driver is illustrated and its command-line outputs are shown. 
```bash
(base) Kathryns-MacBook-Pro-2:example kmaterna$ driver_single_plot.py P511
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
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/drivers_and_configs/example_pngs/P511_noseasons_lssq_cwu_NA_ts.jpg)

### GNSS Stacks
We can also plot a list of GNSS time series as a stack.  An example is shown below for stations in the northern San Francisco Bay Area. The parameters can be set in the 'driver_stack.py' file. 
![GNSS_TS](https://github.com/kmaterna/Mendocino_Geodesy/blob/master/drivers_and_configs/example_pngs/Nbay_-124.0_38.0_125_TS_noeq.jpg)



## Contributing
If you use this library, please contribute back any features you'd like to write! 


## Dependencies and Installation
This library requires basic Python packages such as numpy and matplotlib. For the "stations_within_radius" function, it also requires a file called 'haversine' that is located in my Utilities library (https://github.com/kmaterna/Utility_Code), so please download or clone that as well. 

Please download, clone, or fork this library and put the path to GPS_TOOLS/ on your $PYTHONPATH. The drivers in drivers_and_configs should provide some examples on how to access the major functions in this library.  I expect this library would work on Mac and Linux systems, but possibly not Windows. 


## Acessing Prerequisite GNSS Data
This library requires a local copy of GNSS time series, velocities, steps, and hydrological loading models from various online repositories in their online-provided formats. The "Data Instructions" file tells you how to locate and download these files.  Following that, you must create a file called data_config.txt that tells the code in this library where all of those files are located.  An example data_config.txt is provided; please follow a similar template. The path to the data_config.txt file will be passed into the library each time you use it.  

