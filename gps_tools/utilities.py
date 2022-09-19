

def check_lon_sanity(lon):
    """
    Make sure longitude is between -180 and +180
    """
    if lon > 180:
        lon = lon - 360.0;
    if lon < -360:
        lon = lon + 360;
    return lon;


def remove_blacklist(stations, blacklisted_stations):
    """
    :param stations: list of strings
    :param blacklisted_stations: list of strings
    :return: list of strings
    """
    new_stations = [];
    for station in stations:
        if not (station in blacklisted_stations):
            new_stations.append(station);
        else:
            print("Excluding station %s due to blacklist." % station);
    return new_stations;