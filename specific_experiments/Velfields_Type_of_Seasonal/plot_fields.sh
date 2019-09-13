#!/bin/bash

# Plot that shows MTJ velocities and velocity differences between seasonal removal techniques
# 

lonW="-126.0"
lonE="-121.7"  
latS="38.7"  
latN="42.1"
range="$lonW/$lonE/$latS/$latN"
projection="M6.0i"  # used for medium experiments.
horiz_scale=5.9  # used for velocity vectors
output1='Vel_Differences_Hztl.ps'
output2='Vels_Horizontal.ps'
output4='Vels_Vertical.ps'

folder="Fields/"
file1=$folder"pbo_NA_none_small_velocities.txt"
file2=$folder"pbo_NA_lssq_small_velocities.txt"
file3=$folder"pbo_NA_lsdm_small_velocities.txt"
file4=$folder"pbo_NA_nldas_small_velocities.txt"
file5=$folder"pbo_NA_gldas_small_velocities.txt"
file6=$folder"pbo_NA_grace_small_velocities.txt"

# Pasting two velocity fields together (so we can subtract them later). 
paste $file2 $file1 > $folder"file21.txt"
paste $file3 $file1 > $folder"file31.txt"
paste $file4 $file1 > $folder"file41.txt"
paste $file5 $file1 > $folder"file51.txt"
paste $file6 $file1 > $folder"file61.txt"

gmt makecpt -T-29000/8000/500 -Cgray -Z > blue_topo.cpt



# First plot: Velocity differences
# set up topography
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output1
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output1
gmt pscoast -R$range -J$projection -Lf-125.3/39.15/39.15/50+jt -N1 -N2 -Wthinner,black -Dh -K -O >> $output1 # the title goes here

# Add the plate boundaries
gmt psxy ../../../Misc/Mapping_Resources/transform.gmt -R$range -J$projection -Wthin,red -K -O >> $output1
gmt psxy ../../../Misc/Mapping_Resources/ridge.gmt -R$range -J$projection -Wthin,red -K -O >> $output1
gmt psxy ../../../Misc/Mapping_Resources/trench.gmt -R$range -J$projection -Wthin,red -Sf1.5/0.6+r+t+o1.8 -K -O >> $output1

# Add PBO velocity vectors
awk '{print $1, $2, $3-$13, $4-$14, 0, 0, $9}' $folder"file21.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinnest,black >> $output1
awk '{print $1, $2, $3-$13, $4-$14, 0, 0, $9}' $folder"file31.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblue+pthickest -Wthinnest,blue >> $output1
awk '{print $1, $2, $3-$13, $4-$14, 0, 0, $9}' $folder"file41.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinnest,red >> $output1
awk '{print $1, $2, $3-$13, $4-$14, 0, 0, $9}' $folder"file51.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+ggreen+pthickest -Wthinnest,green >> $output1
awk '{print $1, $2, $3-$13, $4-$14, 0, 0, $9}' $folder"file61.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gdarkorchid2+pthickest -Wthinnest,darkorchid2 >> $output1


gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wthinnest,black -K -O <<EOF >> $output1
-125.0 39.4 .1 0 0.01 0.01 0.0 0.1+-0.01 mm/yr
EOF

gmt pslegend -R$range -J$projection -F+gazure1 -Dx0.1i/6.1i+w2i/1.7i+jTL+l1.2 -C0.1i/0.1i -K -O << EOF >> $output1
# Legend test for pslegend
# G is vertical gap, V is vertical line, N sets # of columns, D draws horizontal line.
# H is header, L is label, S is symbol, T is paragraph text, M is map scale.
#
G -0.1i
H 16 Times-Roman Velocities
D 0.2i 1p
N 1
S 0.1i v0.1i+a40+e 0.2i black 1.8p 0.3i LSSQ-NONE
S 0.1i v0.1i+a40+e 0.2i blue 1.8p,blue 0.3i LSDM-NONE
S 0.1i v0.1i+a40+e 0.2i red 1.8p,red 0.3i NLDAS-NONE
S 0.1i v0.1i+a40+e 0.2i green 1.8p,green 0.3i GLDAS-NONE
S 0.1i v0.1i+a40+e 0.2i darkorchid2 1.8p,darkorchid2 0.3i GRACE-NONE
D 0.2i 1p
P
EOF



# Second plot: The overall velocities
horiz_scale=0.05
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output2
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output2
gmt pscoast -R$range -J$projection -Lf-125.3/39.15/39.15/50+jt -N1 -N2 -Wthinner,black -Dh -K -O >> $output2 # the title goes here

# Add the plate boundaries
gmt psxy ../../../Misc/Mapping_Resources/transform.gmt -R$range -J$projection -Wthin,red -K -O >> $output2
gmt psxy ../../../Misc/Mapping_Resources/ridge.gmt -R$range -J$projection -Wthin,red -K -O >> $output2
gmt psxy ../../../Misc/Mapping_Resources/trench.gmt -R$range -J$projection -Wthin,red -Sf1.5/0.6+r+t+o1.8 -K -O >> $output2

# Add PBO velocity vectors
awk '{print $1, $2, $3, $4, 0, 0, $9}' $folder"file21.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblack+pthickest -Wthinnest,black >> $output2
awk '{print $1, $2, $3, $4, 0, 0, $9}' $folder"file31.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gblue+pthickest -Wthinnest,blue >> $output2
awk '{print $1, $2, $3, $4, 0, 0, $9}' $folder"file41.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gred+pthickest -Wthinnest,red >> $output2
awk '{print $1, $2, $3, $4, 0, 0, $9}' $folder"file51.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+ggreen+pthickest -Wthinnest,green >> $output2
awk '{print $1, $2, $3, $4, 0, 0, $9}' $folder"file61.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+e+gdarkorchid2+pthickest -Wthinnest,darkorchid2 >> $output2

gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wthinnest,black -K -O <<EOF >> $output2
-125.0 39.4 20 0 2 2 0.0 20+-2 mm/yr
EOF

gmt pslegend -R$range -J$projection -F+gazure1 -Dx0.1i/6.1i+w2i/1.7i+jTL+l1.2 -C0.1i/0.1i -K -O << EOF >> $output2
# Legend test for pslegend
# G is vertical gap, V is vertical line, N sets # of columns, D draws horizontal line.
# H is header, L is label, S is symbol, T is paragraph text, M is map scale.
#
G -0.1i
H 16 Times-Roman Velocities
D 0.2i 1p
N 1
S 0.1i v0.1i+a40+e 0.2i black 1.8p 0.3i LSSQ
S 0.1i v0.1i+a40+e 0.2i blue 1.8p,blue 0.3i LSDM
S 0.1i v0.1i+a40+e 0.2i red 1.8p,red 0.3i NLDAS
S 0.1i v0.1i+a40+e 0.2i green 1.8p,green 0.3i GLDAS
S 0.1i v0.1i+a40+e 0.2i darkorchid2 1.8p,darkorchid2 0.3i GRACE
D 0.2i 1p
P
EOF



# The vert velocities with MIDAS
horiz_scale=0.3
gmt makecpt -T-3/3/0.05 -Cpolar -D > vert.cpt
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output4
gmt grdgradient ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -A320 -R$range -Getopo1.grad -Nt
gmt grdhisteq etopo1.grad -Getopo1.hist -N
gmt grdinfo etopo1.hist 
gmt grdmath etopo1.hist 8.41977 DIV = etopo1.norm
gmt grdimage ../../../Misc/Mapping_Resources/Global_topography_data/ETOPO1_Bed_g_gmt4.grd -Ietopo1.norm -R$range -J$projection -Cblue_topo.cpt -K -O >> $output4
gmt pscoast -R$range -J$projection -Lf-125.3/39.15/39.15/50+jt -N1 -N2 -Wthinner,black -Dh -K -O >> $output4 # the title goes here

# Add the plate boundaries
gmt psxy ../../../Misc/Mapping_Resources/transform.gmt -R$range -J$projection -Wthin,red -K -O >> $output4
gmt psxy ../../../Misc/Mapping_Resources/ridge.gmt -R$range -J$projection -Wthin,red -K -O >> $output4
gmt psxy ../../../Misc/Mapping_Resources/trench.gmt -R$range -J$projection -Wthin,red -Sf1.5/0.6+r+t+o1.8 -K -O >> $output4

# Add PBO velocity vectors
awk '{print $1, $2, $5}' $folder"file21.txt" | gmt psxy -R$range -J$projection -Sc0.3 -Cvert.cpt -O -K >> $output4  # the colored dot

# The little blades of grass
awk '{print $1-0.02, $2, 0, $5, 0, 0, $9}' $folder"file21.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+gblack+pthickest -Wthinnest,black >> $output4
awk '{print $1-0.01, $2, 0, $5, 0, 0, $9}' $folder"file31.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+gblue+pthickest -Wthinnest,blue >> $output4
awk '{print $1+0.00, $2, 0, $5, 0, 0, $9}' $folder"file41.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+gred+pthickest -Wthinnest,red >> $output4
awk '{print $1+0.01, $2, 0, $5, 0, 0, $9}' $folder"file51.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+ggreen+pthickest -Wthinnest,green >> $output4
awk '{print $1+0.02, $2, 0, $5, 0, 0, $9}' $folder"file61.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+gdarkorchid2+pthickest -Wthinnest,darkorchid2 >> $output4
awk '{print $1+0.03, $2, 0, $5, 0, 0, 0}' $folder"MIDAS.txt" | gmt psvelo -R$range -J$projection -O -K -Se$horiz_scale/0.68/8 -A+gorange1+pthickest -Wthinnest,orange1 >> $output4


gmt psvelo -R$range -J$projection -A+e+gblack+pthickest -Se$horiz_scale/0.68/10 -Wthinnest,black -K -O <<EOF >> $output4
-125.0 39.4 1 0 .2 .2 0.0 1+-.2 mm/yr
EOF

gmt psscale -R$range -J$projection -D6.5i/3i+w7/0.7 -S -B0.5:"Lssq uplift":/:"mm/yr": -Cvert.cpt -K -O >> $output4
# -B$legend_labelstep:$legend_text:/:$legend_unit:

gmt pslegend -R$range -J$projection -F+gazure1 -Dx0.1i/6.1i+w2i/1.7i+jTL+l1.2 -C0.1i/0.1i -K -O << EOF >> $output4
# Legend test for pslegend
# G is vertical gap, V is vertical line, N sets # of columns, D draws horizontal line.
# H is header, L is label, S is symbol, T is paragraph text, M is map scale.
#
G -0.1i
H 16 Times-Roman Velocities
D 0.2i 1p
N 1
S 0.1i v0.1i+a40+e 0.2i black 1.8p 0.3i LSSQ
S 0.1i v0.1i+a40+e 0.2i blue 1.8p,blue 0.3i LSDM
S 0.1i v0.1i+a40+e 0.2i red 1.8p,red 0.3i NLDAS
S 0.1i v0.1i+a40+e 0.2i green 1.8p,green 0.3i GLDAS
S 0.1i v0.1i+a40+e 0.2i darkorchid2 1.8p,darkorchid2 0.3i GRACE
S 0.1i v0.1i+a40+e 0.2i orange1 1.8p,orange1 0.3i MIDAS
D 0.2i 1p
P
EOF


rm etopo1.*
rm gmt.history

open $output4