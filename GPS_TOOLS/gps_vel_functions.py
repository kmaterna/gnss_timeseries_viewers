# A library of functions that operate on velocity-field objects

import gps_io_functions

# For reference: Velfield = collections.namedtuple("Velfield", ['name', 'nlat', 'elon', 'n', 'e', 'u', 'sn', 'se', 'su',
#                                                   'first_epoch', 'last_epoch']);  # in mm/yr, with -180<lon<180


def clean_velfield(velfield, num_years=0, max_horiz_sigma=1000, max_vert_sigma=1000, coord_box=(-180, 180, -90, 90)):
    # Take the raw GPS velocities, and clean them up.
    # Remove data that's less than num_years long,
    # has formal uncertainties above max_sigma,
    # or is outside our box of interest.
    # Default arguments are meant to be global.
    name = [];
    nlat, elon = [], [];
    n, e, u = [], [], [];
    sn, se, su = [], [], [];
    first_epoch = [];
    last_epoch = [];
    for i in range(len(velfield.n)):
        if velfield.sn[i] > max_horiz_sigma:  # too high sigma, please exclude
            continue;
        if velfield.se[i] > max_horiz_sigma:
            continue;
        if velfield.su[i] > max_vert_sigma:
            continue;
        deltatime = velfield.last_epoch[i] - velfield.first_epoch[i];
        if deltatime.days <= num_years * 365.24:  # too short time interval, please exclude
            continue;
        if coord_box[0] < velfield.elon[i] < coord_box[1] and coord_box[2] < velfield.nlat[i] < coord_box[3]:
            # The station is within the box of interest.
            name.append(velfield.name[i]);
            nlat.append(velfield.nlat[i]);
            elon.append(velfield.elon[i]);
            n.append(velfield.n[i]);
            sn.append(velfield.sn[i]);
            e.append(velfield.e[i]);
            se.append(velfield.se[i]);
            u.append(velfield.u[i]);
            su.append(velfield.su[i]);
            first_epoch.append(velfield.first_epoch[i]);
            last_epoch.append(velfield.last_epoch[i]);
    myVelfield = gps_io_functions.Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su,
                                           first_epoch=first_epoch, last_epoch=last_epoch);
    return [myVelfield];


def remove_duplicates(velfield):
    name = [];
    nlat, elon = [], [];
    n, e, u = [], [], [];
    sn, se, su = [], [], [];
    first_epoch = [];
    last_epoch = [];

    for i in range(len(velfield.n)):
        is_duplicate = 0;
        for j in range(len(name)):
            if abs(nlat[j] - velfield.nlat[i]) < 0.0005 and abs(elon[j] - velfield.elon[i]) < 0.0005:
                # we found a duplicate measurement.
                is_duplicate = 1;
            # Right now assuming all entries at the same lat/lon have the same velocity values.

        if is_duplicate == 0:
            name.append(velfield.name[i]);
            nlat.append(velfield.nlat[i]);
            elon.append(velfield.elon[i]);
            n.append(velfield.n[i]);
            sn.append(velfield.sn[i]);
            e.append(velfield.e[i]);
            se.append(velfield.se[i]);
            u.append(velfield.u[i]);
            su.append(velfield.su[i]);
            first_epoch.append(velfield.first_epoch[i]);
            last_epoch.append(velfield.last_epoch[i]);

    myVelfield = gps_io_functions.Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su,
                                           first_epoch=first_epoch, last_epoch=last_epoch);
    return [myVelfield];


def pair_velocity_fields():
    return;


def combine_several_velocity_fields():
    return;
