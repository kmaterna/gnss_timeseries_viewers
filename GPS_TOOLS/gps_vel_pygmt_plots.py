# Make plots of velocity fields with PyGMT

import numpy as np
import pygmt


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


def simple_pygmt_plot(velfield, outname, symsize=0.1, region=(), horiz_velfield=None):
    """Simply plot the displacement vectors from a velocity field in PyGMT, with vertical denoted by colors"""
    lons = [x.elon for x in velfield];
    lats = [x.nlat for x in velfield];
    erange = np.max(lons) - np.min(lons);
    nrange = np.max(lats) - np.min(lats);
    if not region:
        region = [np.min(lons)-erange*0.06, np.max(lons)+erange*0.06,
                  np.min(lats)-nrange*0.06, np.max(lats)+nrange*0.06];

    fig = pygmt.Figure()
    pygmt.makecpt(cmap="polar", series="-5/5/0.25", background="o", output="mycpt.cpt");
    fig.coast(region=region, projection="M7i", frame="0.5", borders='2', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50");
    elon, nlat, e, n, u = station_vels_to_arrays(velfield);
    fig.plot(x=elon, y=nlat, style='c0.04i', color='black', pen='0.4p,white');  # station locations
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
             direction=[e, n]);  # displacement vectors
    if horiz_velfield:  # horizontal-only velocities
        h_elon, h_nlat, h_e, h_n, h_u = station_vels_to_arrays(horiz_velfield);
        fig.plot(x=h_elon, y=h_nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
                 direction=[h_e, h_n]);  # displacement vectors
    fig.plot(x=elon, y=nlat, style='c'+str(symsize)+'i', G=u, cmap='mycpt.cpt', pen='thin,black');  # vertical
    fig.plot(x=region[0] + 0.9, y=region[2] + 0.1, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
             direction=[[20], [0]]);  # scale vector
    fig.text(x=region[0] + 0.5, y=region[2] + 0.1, text="20 mm/yr", font='10p,Helvetica,black')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", frame=["x1", "y+L\"mm/yr\""]);
    fig.savefig(outname);
    print("Saving %s " % outname);
    return;
