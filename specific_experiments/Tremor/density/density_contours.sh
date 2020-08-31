#!/bin/bash

lonW=-125.5
lonE=-121.2
latS=39.0
latN=42.1

range="$lonW/$lonE/$latS/$latN"
projection="M5.5i"
output="Mend_density.ps"

gmt makecpt -T-2000/2000/200 -Cjet > myrainbow.cpt

# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Lf-121.7/39.25/39.25/50+jt -Wblack -N1 -N2 -B1.0:."Tremor at MTJ": -Dh -K > $output 

gmt grdcontour density.nc -Cmyrainbow.cpt -A- -W+cl -R$range -Wthick -J$projection -O >> $output
# -W+cl takes the colors of the contour lines from the cpt

open $output

gmt psconvert $output -Tg

rm gmt.history
rm myrainbow.cpt