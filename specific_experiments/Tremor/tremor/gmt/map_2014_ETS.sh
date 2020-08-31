#! /bin/bash
# Make an image of the subducting slab in Mendocino
# 10/04/2016 Kathryn Materna

lonW=-125.5
lonE=-121.2
latS=39.0
latN=42.1
range="$lonW/$lonE/$latS/$latN"
projection="m1.5i"
bigrange="-127/-114.0/35/46.5"
bigprojection="m0.09i"
output="mendocino_feb2014.ps"
resampling="0.004"  # smaller means more fine sampling (bigger file)

# Make color scale
gmt makecpt -T0/30/2 -Cjet -Z > myrainbow.cpt

# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Lf-121.7/39.15/39.15/50+jt -Wblack -N2 -B1.0:."Tremor at MTJ": -Dh -K -P > $output 

# Take slab from USGS model
gmt grdcontour cas_slab1.0_clip.grd -R$range -J$projection -C10 -Gd8 -A -Wthick,black -K -P -O >> $output

# Plot tremor epicenters
awk '{print $4, $3}' ../../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt | gmt psxy -R$range -J$projection -Sc0.02 -Gblack -K -O -P >> $output
# Tremor on the day of the March 10 earthquake.
awk '{print $1, $2, $3}' ../Feb2014.txt | gmt psxy -R$range -J$projection -Sc0.10 -Cmyrainbow.cpt -K -O -P >> $output  


# Plot boxes for more in-depth analysis
gmt psxy -R$range -J$projection -Wthicker,cyan4 -K -O -P <<EOF >> $output
-124.0 41.0
-123.35 41.0
-123.35 40.0
-124.0 40.0
-124.0 41.0
EOF
gmt psxy -R$range -J$projection -Wthicker,darkorchid1 -K -O -P <<EOF >> $output
-123.3 41.0
-123.0 41.0
-123.0 40.0
-123.3 40.0
-123.3 41.0
EOF
gmt psxy -R$range -J$projection -Wthicker,darkorange1 -K -O -P <<EOF >> $output
-122.9 41.0
-122.0 41.0
-122.0 40.0
-122.9 40.0
-122.9 41.0
EOF


# # Add PBO velocity vectors
# awk '{print $1, $2, $3, $4, $7, $8}' MTJ_2014_lssq.txt | gmt psvelo -R$range -J$projection -O -P -K -Se0.3/0.68/8 -A+e+gblack+pthickest -Wthick,black >> $output
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se0.3/0.68/10 -Wblack -K -O <<EOF >> $output
# -124.85 41.8 1 0 0.0 0.0 0.0 1mm/yr
# EOF

gmt psmeca -R$range -J$projection -Gdarkmagenta -Sm0.5 -C -K -O <<EOF>> $output
-125.13383 40.82867 15 -0.06 -2.84 2.90 0.21 -0.08 0.48 26 0 0 2014M6.8
EOF


#  Put a tiny map in the corner to orient you
gmt pscoast -R$bigrange -J$bigprojection -Ggray -SWhite -Di -N2 -K -O -P >> $output
gmt psxy -R$bigrange -J$bigprojection -Wthick,red -K -O -P <<EOF >> $output
$lonW $latS
$lonE $latS
$lonE $latN
$lonW $latN
$lonW $latS
EOF

# Plot the scale for the slip rate colors
gmt psscale -R -J -DjTR+w3.4c/0.5c+o-1.5/-0.3 -Cmyrainbow.cpt -B10.0:"Days":/:: -P -O -K >> $output
# D: Dimensions (MUST BE IN inches / centimeters)
# B: Scale has 5km boxes.  I just copied the syntax on the annotation.  

