# A library of functions that operate on velocity-field objects (lists of station-vels)


def clean_velfield(velfield, num_years=0, max_horiz_sigma=1000, max_vert_sigma=1000, coord_box=(-180, 180, -90, 90)):
    # Filter into cleaner GPS velocities:
    # Remove data that's less than num_years long, has formal uncertainties above max_sigma,
    # or is outside our box of interest.
    # Default arguments are meant to be global.
    cleaned_velfield = [];
    for station_vel in velfield:
        if station_vel.sn > max_horiz_sigma:  # too high sigma, please exclude
            continue;
        if station_vel.se > max_horiz_sigma:
            continue;
        if station_vel.su > max_vert_sigma:
            continue;
        deltatime = station_vel.last_epoch - station_vel.first_epoch;
        if deltatime.days <= num_years * 365.24:  # too short time interval, please exclude
            continue;
        if coord_box[0] < station_vel.elon < coord_box[1] and coord_box[2] < station_vel.nlat < coord_box[3]:
            # The station is within the box of interest.
            cleaned_velfield.append(station_vel);
    return cleaned_velfield;


def remove_duplicates(velfield):
    # Right now assuming all entries at the same lat/lon are the same station
    cleaned_velfield = []
    for vel in velfield:
        is_duplicate = 0;
        for comp_station_vel in cleaned_velfield:
            if abs(comp_station_vel.elon - vel.elon) < 0.0005 and abs(comp_station_vel.nlat - vel.nlat) < 0.0005:
                # we found a duplicate measurement.
                is_duplicate = 1;
        if is_duplicate == 0:
            cleaned_velfield.append(vel);
    return cleaned_velfield;


def remove_blacklist_vels(velfield, blacklist):
    cleaned_velfield = []
    for vel in velfield:
        if vel.name not in blacklist:
            cleaned_velfield.append(vel);
    return cleaned_velfield;


def pair_velocity_fields():
    return;


def combine_several_velocity_fields():
    return;
