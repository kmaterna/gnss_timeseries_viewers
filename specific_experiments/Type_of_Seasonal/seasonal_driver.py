#   Goal: Develop tools for removing seasonals by...
#    fit: fits seasonals and linear trend by least squares inversion.
#   noel: uses noel's fits of inter-SSE velocities and seasonal terms.
#  notch: removes the 1-year and 6-month components by notch filter.
#  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not available.
#    stl: not supported yet. 

# Goal 1: Implement these seasonal corrections for a single station. Plot them all next to each other. Done. 
# Goal 2: Implement these for a bunch of stations, map the timing. Done. 
# Goal 3: Plot the stacked time series with a "corner" model on top, once for each seasonal method. Optional. 


import seasonal_single_plot


station="P162";
seasonal_single_plot.compare_single_seasonals(station, offsets_remove=1, earthquakes_remove=1, outliers_remove=1);

