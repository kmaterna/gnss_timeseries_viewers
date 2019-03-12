#! /bin/bash
# Plot GPS Vectors and Topography and Faults at Mendocino Triple Junction. 

infile="OUTPUTS/model_vs_data.txt"
lonW=-125
lonE=-121
latS=39.3
latN=42.1
output1="OUTPUTS/loading_horiz.ps"
horiz_scale=0.25  # used for velocity change vectors (0.3 sometimes, sometimes smaller)
range="$lonW/$lonE/$latS/$latN"
projection="M4.0i"  # used for medium experiments.

gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt

# The Horizontals
# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0wESn:."2014": -Dh -K > $output1 # the title goes here

gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-124.6/39.50/39.50/50+jt -N1 -N2 -Wthinner,black -Dh -K -O >> $output1 # the title goes here

# Add velocity vectors
# KEY: lon, lat, T3-T2_gps(E,N) T3-T2_model(E,N) T4-T3_gps(E,N) T4-T3_model(E,N)
awk '{print $1, $2, $3, $4, 0, 0, 0}' $infile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthick,black >> $output1
awk '{print $1, $2, $5, $6, 0, 0, 0}' $infile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthick,red >> $output1

gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.0 39.65 2 0 0 0 0.0 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-124.0 39.8 2 0 0 0 0.0 2mm/yr model
EOF


# The Horizontals
# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0wESn:."2016": -Dh -K -O -X13 >> $output1 # the title goes here

gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-124.6/39.50/39.50/50+jt -N1 -N2 -Wthinner,black -Dh -K -O >> $output1 # the title goes here

# Add velocity vectors
# KEY: lon, lat, T3-T2_gps(E,N) T3-T2_model(E,N) T4-T3_gps(E,N) T4-T3_model(E,N)
awk '{print $1, $2, $7, $8, 0, 0, 0}' $infile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthick,black >> $output1
awk '{print $1, $2, $9, $10, 0, 0, 0}' $infile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthick,red >> $output1

gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.0 39.65 2 0 0 0 0.0 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-124.0 39.8 2 0 0 0 0.0 2mm/yr model
EOF

rm gmt.history
rm etopo1.grad
rm etopo1.hist
rm etopo1.norm
rm blue_topo.cpt

open $output1