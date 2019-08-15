# This is the plot for the Figure in the paper. 

# The preferred method for the paper: 
# infolder="HET/pbo_nldas_NA_characteristic/"

infolder="HET/pbo_nldas_NA_characteristic/"
name1="2014"
name2="2016"
infile1=$infolder$name1"_inversion.txt"
infile2=$infolder$name2"_inversion.txt"
gpsfile1=$infolder$name1".txt"
gpsfile2=$infolder$name2".txt"
modelfile1=$infolder$name1"_model.txt"
modelfile2=$infolder$name2"_model.txt"
lonW=-125
lonE=-121
latS=39
latN=42.2
output1=$infolder"output.ps"
horiz_scale=0.25  # used for velocity change vectors (0.3 sometimes, sometimes smaller)
# horiz_scale=0.9  # used for velocity change vectors (used for tiny grace arrows)
range="$lonW/$lonE/$latS/$latN"
projection="M4.5i"  # used for medium experiments.

# Make colorscales
# gmt makecpt -T0/29/1 -Ic -Chot > datacpt.cpt
gmt makecpt -T-22/22/1 -CBlueWhiteOrangeRed > datacpt2.cpt # this is for the 2016 case
gmt makecpt -T-22/22/1 -Ic -CBlueWhiteOrangeRed > datacpt2_backwards.cpt  # this is for the 2014 case
gmt makecpt -T-22/22/1 -CBlueWhiteOrangeRed > datacpt2_scale.cpt # this is for the 2014 color scale
gmt makecpt -T-1000/2000/300 -Ic -Ccopper > tremor.cpt  # Cdrywet also works
gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt


# gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile1: -Dh -K > $output1 # the title goes here
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K > $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.19/39.19/50+jt -Spaleturquoise -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile1 -R$range -J$projection -L -Wthinner,gray -Cdatacpt2_backwards.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1


gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1



# -W+cl takes the colors of the contour lines from the cpt
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile1 | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthin,black >> $output1
awk '{print $1, $2, $3*-1, $4*-1}' $modelfile1 | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthin,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.1 39.4 2 0 0.0 0.0 0.0 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-124.1 39.2 2 0 0.0 0.0 0.0 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
-121.7 42.1 $name1: T3-T2
-123.9 42.1 A: Coupling Increase
EOF
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se0.04/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.9 41.75 28 14 0 0 0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
-124.58 41.7 31 mm/yr
-124.58 41.6 plate
-124.55 41.5 convergence
EOF


# The second image
# gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0weSn:.$infile2: -Dh -K -O -X12.5 >> $output1 # the title goes here
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0weSn -Dh -K -O -X12.5 >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.19/39.19/50+jt -Spaleturquoise -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile2 -R$range -J$projection -L -Wthinner,gray -Cdatacpt2.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1

awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile2 | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthin,black >> $output1
gmt psvelo $modelfile2 -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthin,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-124.1 39.4 2 0 0.0 0.0 0.0 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-124.1 39.2 2 0 0.0 0.0 0.0 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
-121.7 42.1 $name2: T4-T3
-123.9 42.1 B: Coupling Decrease
EOF

# gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-1.5/0.2 -Cdatacpt.cpt -B5.0:"":/:mm/yr: -P -O -K >> $output1
gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-1.5/0.2 -G-22/22 -Cdatacpt2.cpt -B5.0:"":/:mm/yr: -P -O >> $output1
# gmt psscale -R$range -J$projection -DjTR+w6c/0.5c+o-1.5/6.2 -G-22/-1 -Cdatacpt2_scale.cpt -B5.0:"":/:"": -P -O >> $output1

gmt psconvert $output1 -Tg

open $output1

