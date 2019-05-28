#!/bin/bash

lonE=-119.8
lonW=-125
latS=39
latN=42.1
projection=M5i
range=$lonW/$lonE/$latS/$latN
output="box.ps"

gmt makecpt -T0/40/1 -Cjet > mycpt.cpt

gmt pscoast -R$range -J$projection -Wthin,black -Swhite -Gwhite -B1.0 -N1 -N2 -Dh -P -K > $output

awk '{print $4, $3, $5, $6*0.05}' ../../GPS_POS_DATA/Seismicity/MTJ_box_search.txt | gmt psxy -R$range -J$projection -Sc -Cmycpt.cpt -P -K -O >> $output

gmt psscale -R$range -J$projection -DjTR+w6c/0.5c+o-1.5/0.2 -Cmycpt.cpt -B5.0:"":/:km: -P -K -O >> $output

gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthin,purple -C30 -A -K -O >> $output

gmt psxy -R$range -J$projection -Wthick,red -P -K -O <<EOF>> $output
-123.6 40.8 
-122.8 40.8
-122.8 40.0
-123.6 40.0
-123.6 40.8
EOF

open $output