import collections


# Station_Vel are in mm/yr, with -180<lon<180, used for velfields
Station_Vel = collections.namedtuple("Station_Vel", ['name', 'nlat', 'elon', 'n', 'e', 'u', 'sn', 'se', 'su',
                                                     'first_epoch',
                                                     'last_epoch', 'refframe', 'proccenter', 'subnetwork', 'survey',
                                                     'meas_type'], defaults=(None,) * 15 + ('gnss',));

# Station_Vel_XYZ are in m/yr (under consideration), with XYZ as ECEF position in meters. Used for velocities.
Station_Vel_XYZ = collections.namedtuple("Station_Vel_XYZ", ['name',
                                                             'x_pos', 'y_pos', 'z_pos',
                                                             'x_rate', 'y_rate', 'z_rate',
                                                             'x_sigma', 'y_sigma', 'z_sigma',
                                                             'first_epoch', 'last_epoch']);

Timeseries = collections.namedtuple("Timeseries", ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su',
                                                   'EQtimes']);  # in mm


# The namedtuple definition.  Offsets should be in mm. One object per offset. Each field is a single value.
Offset = collections.namedtuple("Offsets", ['e_offset', 'n_offset', 'u_offset', 'evdt']);


# u, v, w are GRACE model displacements in east, north, and up.
Paired_TS = collections.namedtuple('Paired_TS', ['dtarray', 'north', 'east', 'vert', 'N_err', 'E_err', 'V_err',
                                                 'u', 'v', 'w']);
