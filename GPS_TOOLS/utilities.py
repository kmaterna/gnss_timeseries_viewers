

def check_lon_sanity(lon):
    """
    Make sure longitude is between -180 and +180
    """
    if lon > 180:
        lon = lon - 360.0;
    if lon < -360:
        lon = lon + 360;
    return lon;
