"""Make plots of velocity fields with PyGMT"""

import numpy as np
import pygmt
from . import vel_functions


def station_vels_to_arrays(station_vels):
    """ Unpack station vels into arrays for pygmt plotting of vectors """
    elon, nlat, e, n, u = [0], [0], [0], [0], [0];
    for item in station_vels:
        elon.append(item.elon);
        nlat.append(item.nlat);
        e.append(item.e)
        n.append(item.n);
        u.append(item.u);
    return np.array(elon), np.array(nlat), np.array(e), np.array(n), np.array(u);


def simple_pygmt_plot(velfield, outname, symsize=0.1, region=(), horiz_velfield=None, fault_file=None,
                      vector_scale_info=(1, 2, "2 mm"), plot_names=True):
    """
    Simply plot the displacement vectors from StationVels in PyGMT, with vertical denoted by colors
    example: scale_factor = 1; scale_arrow = 2; vectext = "2 mm";  # 1 mm, extra small vectors
    """
    lons = [x.elon for x in velfield];
    lats = [x.nlat for x in velfield];
    names = [x.name for x in velfield];
    erange = np.max(lons) - np.min(lons);
    nrange = np.max(lats) - np.min(lats);
    if not region:
        region = [np.min(lons)-erange*0.06, np.max(lons)+erange*0.06,
                  np.min(lats)-nrange*0.06, np.max(lats)+nrange*0.06];

    fig = pygmt.Figure()
    pygmt.makecpt(cmap="polar", series="-5/5/0.25", background="o", output="mycpt.cpt");
    fig.coast(region=region, projection="M7i", frame="0.5", borders=['1', '2'], shorelines='1.0p,black',
              water='lightblue', map_scale="n0.12/0.12+c" + str(region[2]) + "+w10");
    if fault_file:
        fig.plot(data=fault_file, pen='thin,darkgray');
    elon, nlat, e, n, u = station_vels_to_arrays(velfield);
    scale = (6 * vector_scale_info[0] / 12);  # empirical scaling, convenient display
    fig.plot(x=elon, y=nlat, style='c0.04i', color='black', pen='0.4p,white');  # station locations
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z'+str(scale), pen='0.6p,black',
             direction=[e, n]);  # displacement vectors
    if horiz_velfield:  # horizontal-only velocities
        h_elon, h_nlat, h_e, h_n, h_u = station_vels_to_arrays(horiz_velfield);
        fig.plot(x=h_elon, y=h_nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z'+str(scale), pen='0.6p,black',
                 direction=[h_e, h_n]);  # displacement vectors
    fig.plot(x=elon, y=nlat, style='c'+str(symsize)+'i', color=u, cmap='mycpt.cpt', pen='thin,black');  # vertical
    fig.plot(x=region[0], y=region[2], style='v0.20+e+a40+gblack+h0+p1p,black+z'+str(scale), offset="0.5i/0.3i",
             pen='0.6p,black', direction=[[vector_scale_info[1]], [0]]);  # scale vector
    fig.text(x=region[0], y=region[2], text=vector_scale_info[2], offset="0.3i/0.3i", font='10p,Helvetica,black')
    if plot_names:
        fig.text(x=[i + 0.02 for i in lons], y=lats, text=names, font='11p,Helvetica-Bold,black');
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", frame=["x1", "y+L\"mm/yr\""]);
    fig.savefig(outname);
    print("Saving %s " % outname);
    return;


def map_velocity_profile(velfield, selected_velfield, outname,  vector_scale_info=(0.1, 2, "20 mm"),
                         startcoord=None, endcoord=None, fault_traces=None):
    """
    Map a velocity field.
    Then map a selected profile of stations across it, in a different color.
    fault_traces can be a list of files
    """
    fig = pygmt.Figure();
    region = vel_functions.get_bounding_box(velfield);
    fig.coast(region=region, projection="M6i", frame="1.0", shorelines="1.0p,black", water="lightblue",
              borders=['1', '2']);
    elon, nlat, e, n, u = station_vels_to_arrays(velfield);
    scale = (6 * vector_scale_info[0] / 12);  # empirical scaling, convenient display
    fig.plot(x=elon, y=nlat, style='c0.04i', color='black', pen='0.4p,white');  # station locations
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z'+str(scale), pen='0.6p,black',
             direction=[e, n]);  # displacement vectors

    # The selected stations in the profile
    elon, nlat, e, n, u = station_vels_to_arrays(selected_velfield);
    scale = (6 * vector_scale_info[0] / 12);  # empirical scaling, convenient display
    fig.plot(x=elon, y=nlat, style='c0.04i', color='red', pen='0.4p,white');  # station locations
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gred+h0+p1p,red+z'+str(scale), pen='0.6p,red',
             direction=[e, n]);  # displacement vectors

    if fault_traces:
        for item in fault_traces:
            fig.plot(data=item, pen='1.6p,slateblue3');

    if startcoord and endcoord:
        fig.plot(x=[startcoord[0], endcoord[0]], y=[startcoord[1], endcoord[1]], pen='0.4p,black,dashed');

    fig.savefig(outname);
    return;


def map_ts_objects(dataobj_list, outname, center=None):
    """
    Plot stations locations, maybe with a star for the center of the search radius

    :param dataobj_list: list of TS objects
    :param outname: string
    :param center: [lon, lat], optional
    """
    offset = 0.2;

    lons = [x.coords[0] for x in dataobj_list];
    lats = [x.coords[1] for x in dataobj_list];
    names = [x.name for x in dataobj_list];
    region = [min(lons) - offset, max(lons) + offset, min(lats) - offset, max(lats) + offset];

    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M8i", frame="0.25");
    # fig.grdimage("@earth_relief_30s",region=region,I="+d");  # takes a while the first time, but faster afterwards
    fig.coast(shorelines="0.5p,black", land='peachpuff2', water='skyblue', resolution="h");
    fig.coast(borders='1', shorelines='1.0p,black');
    fig.coast(borders='2', shorelines='0.5p,black');
    fig.text(x=lons, y=lats, text=names, font='15p,Helvetica-Bold,black', offset="0.24i/0.11i");
    fig.plot(x=lons, y=lats, style='c0.1i', fill='black', pen='0.5p,black')
    if center:
        fig.plot(x=center[0], y=center[1], style='a0.1i', fill='red', pen='0.5p,red')
    fig.savefig(outname);
    print("Saving map %s" % outname );
    return;
