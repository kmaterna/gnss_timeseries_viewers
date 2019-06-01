# Remove data during ETS times
# Estimate offsets during those times. 

import remove_ets_events 


if __name__=="__main__":
	station="P330"
	remove_ets_events.view_single_station(station, 
		offsets_remove=1, earthquakes_remove=1, 
		outliers_remove=1, seasonals_remove=1,seasonals_type='lssq', datasource='pbo', refframe='NA');



