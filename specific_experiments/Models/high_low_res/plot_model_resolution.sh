# This is the plot for the Figure in the paper. 

infile1="2014_over.txt"
infile2="2014_inversion.txt"
infile3="2014_under.txt"
infile4="2016_over.txt"
infile5="2016_inversion.txt"
infile6="2016_under.txt"


gpsfile="2014.txt"
lonW=-125.0
lonE=-121
latS=39
latN=42.42
output1="resolution.ps"
range="$lonW/$lonE/$latS/$latN"
horiz_scale=0.25  # used for velocity change vectors (0.3 sometimes, sometimes smaller)
projection="M3.0i"  # used for medium experiments.
label_lat="42.3"
watercolor=-Spaleturquoise

# Make colorscales
gmt makecpt -T-26/26/1 -Ic -CBlueWhiteOrangeRed > datacpt.cpt  # this is for the 2014 case
gmt makecpt -T-26/26/1 -CBlueWhiteOrangeRed > datacpt_upside.cpt  # this is for the 2016 case
gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt


# The first plot
# gpsfile=$infolder1$name1".txt"
# modelfile=$infolder1$name1"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X1 -Y12 -K > $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile1 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.15 $label_lat A: 2014 over-smoothed (3200)
EOF





# The second plot
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X9 -K -O >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile2 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.65 $label_lat B: 2014 selected (800)
EOF



# The third plot
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X9 -K -O >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile3 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.15 $label_lat C: 2014 under-smoothed (200)
EOF







# The fourth plot
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X-18 -Y-10 -K -O >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile4 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt_upside.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.15 $label_lat D: 2016 over-smoothed (3200)
EOF


# The fifth plot
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X9 -K -O >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile5 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt_upside.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.65 $label_lat E: 2016 selected (800)
EOF






# The sixth plot
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:."": -Dh -X9 -K -O >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile6 -R$range -J$projection -L -Wthinnest,gray -Cdatacpt_upside.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

# gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
gmt grdcontour ../../Tremor/slab_geometry/cas_slab1.0_clip.grd -R$range -J$projection -Wthicker,purple -C30 -A+f10+ukm+an+ggray80 -GL-124.5/41.8/-120/41.2 -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
# awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
# awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
# gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
# -125.4 39.6 2 0 0.0 0.0 0.0 
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
# -124.9 39.8 2mm/yr data
# EOF
# gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
# -125.4 39.2 2 0 0.0 0.0 0.0
# EOF
# gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
# -124.85 39.4 2mm/yr model
# EOF
gmt pstext -R$range -J$projection -F+f14p,Helvetica -Gwhite -K -O <<EOF >> $output1
-123.15 $label_lat F: 2016 under-smoothed (200)
EOF



gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-2.0/-6.2 -G-26/26 -Cdatacpt.cpt -B5:"":/:mm: -P -O >> $output1



gmt psconvert $output1 -Tg
open $output1

