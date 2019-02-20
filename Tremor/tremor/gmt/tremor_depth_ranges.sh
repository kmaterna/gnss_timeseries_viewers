#! /bin/bash
# Make an image of the subducting slab in Mendocino, coloring tremor by depth boxes

lonW=-125.5
lonE=-121.2
latS=39.0
latN=42.1
range="$lonW/$lonE/$latS/$latN"
projection="m1.5i"
bigrange="-127/-114.0/35/46.5"
bigprojection="m0.09i"
output="MTJ_Tremor_key.ps"

# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Lf-121.7/39.15/39.15/50+jt -Wthick,black -N2 -B1.0:."Tremor at MTJ": -Dh -K -P > $output 

# Plot tremor epicenters
awk '{print $4, $3}' ../../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt | gmt psxy -R$range -J$projection -Sc0.02 -Gblack -K -O -P >> $output

# Plot tremor categories for more in-depth analysis
gmt psxy shallowrange.txt -R$range -J$projection -Sc0.1 -Gcyan4 -K -O -P >> $output
gmt psxy medrange.txt -R$range -J$projection -Sc0.1 -Gdarkorchid1 -K -O -P >> $output
gmt psxy deeprange.txt -R$range -J$projection -Sc0.05 -Gdarkorange1 -K -O -P >> $output

# Slab from USGS model
gmt grdcontour cas_slab1.0_clip.grd -R$range -J$projection -C10 -Gd8 -A -Wthick,black -K -P -O >> $output

# Add PBO velocity vectors
awk '{print $1, $2, $3, $4, $7, $8}' MTJ_2014_lssq.txt | gmt psvelo -R$range -J$projection -O -P -K -Se0.3/0.68/8 -A+e+gblack+pthickest -Wthick,black >> $output
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se0.3/0.68/10 -Wblack -K -O <<EOF >> $output
-124.85 41.8 1 0 0.0 0.0 0.0 1mm/yr
EOF

gmt psmeca -R$range -J$projection -Gblack -Sm0.5 -C -K -O <<EOF>> $output
-125.13383 40.82867 15 -0.06 -2.84 2.90 0.21 -0.08 0.48 26 0 0 2014M6.8
EOF

#  Put a tiny map in the corner to orient you
gmt pscoast -R$bigrange -J$bigprojection -Ggray -SWhite -Di -N2 -K -O -P >> $output
gmt psxy -R$bigrange -J$bigprojection -Wthick,red -O -P <<EOF >> $output
$lonW $latS
$lonE $latS
$lonE $latN
$lonW $latN
$lonW $latS
EOF

gmt psconvert $output -Tg

open $output

