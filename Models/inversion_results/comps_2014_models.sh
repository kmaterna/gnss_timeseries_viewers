# Make 4 models: 
# Heterogeneous Green's functions, ETS corrected
# Homogeneous Green's functions, ETS corrected
# Heterogeneous Green's functions, no ETS correction
# Homogeneous Green's functions, no ETS correction


name1="2014"
name2="2016"
infolder1="HET/pbo_nldas_NA_characteristic/"
infolder2="HET/pbo_nldas_NA_characteristic_homog/"
infolder3="pbo_nldas_NA_HET/"
infolder4="pbo_nldas_NA/"


lonW=-125.7
lonE=-121
latS=39
latN=42.5
output1="GF_comparison_2014.ps"
output2="GF_comparison_2016.ps"
horiz_scale=0.25  # used for velocity change vectors (0.3 sometimes, sometimes smaller)
range="$lonW/$lonE/$latS/$latN"
projection="M3.0i"  # used for medium experiments.
label_lat="42.3"

# Make colorscales
gmt makecpt -T-24/24/1 -Ic -CBlueWhiteOrangeRed > datacpt.cpt  # this is for the 2014 case
gmt makecpt -T-24/24/1 -CBlueWhiteOrangeRed > upside.cpt  # this is for the 2014 case
gmt makecpt -T-40/40/1 -CBlueWhiteOrangeRed > scalecpt.cpt  # this is for the 2016 case

gmt makecpt -T-4000/4000/300 -Ic -Ccopper > tremor.cpt
gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt

watercolor=-Spaleturquoise

# The first plot
infile=$infolder1$name1"_inversion.txt"
gpsfile=$infolder1$name1".txt"
modelfile=$infolder1$name1"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile1: -Dh -X1 -Y12 -K > $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
# -W+cl takes the colors of the contour lines from the cpt
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
# -121.55 $label_lat $name1
-124.0 $label_lat A: HET, ETS Corr
EOF


# The second plot
infile=$infolder2$name1"_inversion.txt"
gpsfile=$infolder2$name1".txt"
modelfile=$infolder2$name1"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X9 >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
# -121.55 $label_lat $name1
-123.6 $label_lat B: HOMOG, ETS Corr
EOF


# The third plot
infile=$infolder3$name1"_inversion.txt"
gpsfile=$infolder3$name1".txt"
modelfile=$infolder3$name1"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X-9 -Y-10 >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
# -121.55 $label_lat $name1
-123.7 $label_lat C: HET, No ETS Corr
EOF


# The fourth plot
infile=$infolder4$name1"_inversion.txt"
gpsfile=$infolder4$name1".txt"
modelfile=$infolder4$name1"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X9 -Y0 >> $output1 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output1 
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cdatacpt.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output1

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output1
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output1
awk '{print $1, $2, $3*-1, $4*-1}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output1
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output1
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output1
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output1
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output1
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output1
# -121.55 $label_lat $name1
-123.4 $label_lat D: HOMOG, No ETS Corr
EOF

gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-1.5/-6 -Cupside.cpt -B5.0:"":/:mm/yr: -P -O >> $output1

gmt psconvert $output1 -Tg






# Making a 2016 plot. 

# The first plot
infile=$infolder1$name2"_inversion.txt"
gpsfile=$infolder1$name2".txt"
modelfile=$infolder1$name2"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile1: -Dh -X1 -Y12 -K > $output2 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output2
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cscalecpt.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output2

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output2
# -W+cl takes the colors of the contour lines from the cpt
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output2
awk '{print $1, $2, $3, $4}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output2
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output2
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output2
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output2
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output2
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output2
# -121.55 $label_lat $name2
-124.0 $label_lat E: HET, ETS Corr
EOF


# The second plot
infile=$infolder2$name2"_inversion.txt"
gpsfile=$infolder2$name2".txt"
modelfile=$infolder2$name2"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X9 >> $output2 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output2
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cscalecpt.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output2

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output2
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output2
awk '{print $1, $2, $3, $4}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output2
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output2
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output2
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output2
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output2
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output2
# -121.55 $label_lat $name2
-123.6 $label_lat F: HOMOG, ETS Corr
EOF


# The third plot
infile=$infolder3$name2"_inversion.txt"
gpsfile=$infolder3$name2".txt"
modelfile=$infolder3$name2"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X-9 -Y-10 >> $output2 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output2
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cscalecpt.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output2

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output2
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output2
awk '{print $1, $2, $3, $4}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output2
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output2
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output2
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output2
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output2
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output2
# -121.55 $label_lat $name2
-123.7 $label_lat G: HET, No ETS Corr
EOF


# The fourth plot
infile=$infolder4$name2"_inversion.txt"
gpsfile=$infolder4$name2".txt"
modelfile=$infolder4$name2"_model.txt"
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn:.$infile2: -Dh -K -O -X9 -Y0 >> $output2 # the title goes here
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Lf-121.6/39.3/39.3/50+jt $watercolor -N1 -N2 -Dh -O -K >> $output2
gmt psxy $infile -R$range -J$projection -L -Wthinnest,gray -Cscalecpt.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Wthicker,black -N1 -N2 -Dh -K -O >> $output2

gmt grdcontour mapping_data/tremor_density.nc -Ctremor.cpt -A- -W+cl -R$range -Wthick -J$projection -K -O >> $output2
awk '{print $1, $2, $3, $4, $7, $8, $10}' $gpsfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinner,black >> $output2
awk '{print $1, $2, $3, $4}' $modelfile | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinner,red >> $output2
gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wblack -K -O <<EOF >> $output2
-125.4 39.6 2 0 0.0 0.0 0.0 
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica -K -O <<EOF >> $output2
-124.9 39.8 2mm/yr data
EOF
gmt psvelo -R$range -J$projection -A+e+gred+pthickest -Se$horiz_scale/0.68/10 -Wred -K -O <<EOF >> $output2
-125.4 39.2 2 0 0.0 0.0 0.0
EOF
gmt pstext -R$range -J$projection -F+f10p,Helvetica,red -K -O <<EOF >> $output2
-124.85 39.4 2mm/yr model
EOF
gmt pstext -R$range -J$projection -F+f18p,Helvetica -Gwhite -K -O <<EOF >> $output2
# -121.55 $label_lat $name2
-123.4 $label_lat H: HOMOG, No ETS Corr
EOF

gmt psscale -R$range -J$projection -DjTR+w12c/0.5c+o-1.5/-6 -Cscalecpt.cpt -B5.0:"":/:mm/yr: -P -O >> $output2

gmt psconvert $output2 -Tg

open $output1

