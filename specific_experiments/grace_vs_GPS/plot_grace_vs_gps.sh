#!/bin/bash

lonW=-125.8
lonE=-117.7
latS=36.8
latN=44.2

datasource="pbo"
hydro_type="grace"
range="$lonW/$lonE/$latS/$latN"
projection="M6.0i"  # used for medium experiments.
ifile=$hydro_type"_vs_gps_amps_"$datasource".txt"
output1=$hydro_type"_vs_gps_vert"$datasource".ps"
output2=$hydro_type"_vs_gps_east"$datasource".ps"

legend_text='AnnualAmp(Vert)'
legend_unit='mm'
legend_min=0
legend_max=9
legend_labelstep=2

gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt
gmt makecpt -T$legend_min/$legend_max/0.1/1 -Cjet > amptopo.cpt

# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0wESn -Dh -K -P > $output1 # the title goes here

gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1

gmt pscoast -R$range -J$projection -Lf-124.8/38.3/38.3/50+jt -N1 -N2 -Wthinner,black -Dh -K -O -P >> $output1 # the title goes here
awk '{print $1, $2, $5}' $ifile | gmt psxy -R$range -J$projection -Sc0.5 -Camptopo.cpt -Wthinner,black -K -O -P >> $output1
awk '{print $1, $2, $8}' $ifile | gmt psxy -R$range -J$projection -Sh0.25 -Camptopo.cpt -Wthinner,black -K -O -P >> $output1
gmt psscale -R$range -J$projection -DjTR+w4c/0.5c+o-1.5/1.5 -B$legend_labelstep:$legend_text:/:$legend_unit: -G$legend_min/$legend_max -Camptopo.cpt -O -K -P >> $output1

# gmt pslegend -R$range -J$projection -O -K -P <<EOF >> $output1
# This is my wishlist
# EOF


legend_text='AnnualAmp(East)'
legend_unit='mm'
legend_min=0
legend_max=2
legend_labelstep=0.5

gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt
gmt makecpt -T$legend_min/$legend_max/0.1/1 -Cjet > amptopo.cpt

# # Insert a horizontal scale bar and title
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0wESn -Dh -K -P > $output2 # the title goes here

gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2

gmt pscoast -R$range -J$projection -Lf-124.8/38.3/38.3/50+jt -N1 -N2 -Wthinner,black -Dh -K -O -P >> $output2
awk '{print $1, $2, $3}' $ifile | gmt psxy -R$range -J$projection -Sc0.5 -Camptopo.cpt -Wthinner,black -K -O -P >> $output2
awk '{print $1, $2, $6}' $ifile | gmt psxy -R$range -J$projection -Sh0.25 -Camptopo.cpt -Wthinner,black -K -O -P >> $output2
gmt psscale -R$range -J$projection -DjTR+w4c/0.5c+o-1.5/1.5 -B$legend_labelstep:$legend_text:/:$legend_unit: -G$legend_min/$legend_max -Camptopo.cpt -O -K -P >> $output2


open $output1
rm gmt.history
rm blue_topo.cpt amptopo.cpt 
rm etopo1.grad etopo1.hist etopo1.norm
