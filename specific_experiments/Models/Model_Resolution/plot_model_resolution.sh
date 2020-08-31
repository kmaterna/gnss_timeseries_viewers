# This is the plot for the Figure in the paper. 

infile1="resolution.txt"
gpsfile="2014.txt"
lonW=-125
lonE=-121
latS=39
latN=42.2
output1="resolution.ps"
range="$lonW/$lonE/$latS/$latN"
horiz_scale=0.25  # used for velocity change vectors (0.3 sometimes, sometimes smaller)
projection="M4.5i"  # used for medium experiments.

# Make colorscales
gmt makecpt -T-0.5/1.2/0.05 -D -Cocean > datacpt2.cpt
# gmt makecpt -T-10/15/1 -Cocean > datacpt2.cpt
gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt

gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."Median Kernel": -Dh -K > $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.19/39.19/50+jt -Spaleturquoise -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile1 -R$range -J$projection -L -Wthinner,gray -Cdatacpt2.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1

awk '{print $1, $2}' $gpsfile | gmt psxy -R$range -J$projection -O -K -Sc0.1 -Gblack -Wthin,black >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.1 39.4 2 0 0.0 0.0 0.0 2mm/yr data
EOF

gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-1.5/0.2 -G0/1.2 -Cdatacpt2.cpt -B0.5:"":/:mm: -P -O >> $output1

gmt psconvert $output1 -Tg
open $output1

