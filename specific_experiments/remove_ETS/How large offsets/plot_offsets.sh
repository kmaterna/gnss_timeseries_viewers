#!/bin/bash

# Plot that shows MTJ data in fancy way

lonW="-125.6"
lonE="-120.3"   # MTJ: -121.2.  NorCal: -120.2
latS="38.2"  # MTJ: 38.7.  NorCal: 36.7
latN="43.5"
range="$lonW/$lonE/$latS/$latN"
projection="M4.5i"  # used for medium experiments.
output1='Offsets.ps'
output2='Offsets_vert.ps'
output3='Offsets_std.ps'
infile=Offsets_2mm_30days.txt
offset=13

gmt makecpt -T-2/2/0.2 -Cjet > colors.cpt

# THE FIRST MAP
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output1
gmt pscoast -R$range -J$projection -Lf-124.9/39.15/39.15/50+jt -N1 -N2 -B+t"Average East Offsets" -Wthinner,black -Dh -K -O >> $output1 # the title goes here

# Add data
awk '{print $2, $3, $4}' $infile | gmt psxy -R$range -J$projection -Sc0.3 -Ccolors.cpt -Wthinnest,black -K -O >> $output1
awk '{print $2, $3, $1}' $infile | gmt pstext -R$range -J$projection -D0/.15 -F+f6p,Helvetica -K -O >> $output1
awk '{print $2, $3, $4, $5}' $infile | gmt psvelo -R$range -J$projection -Se0.2/0.68/0 -A+e+gblack+pthickest -Wthick,black -K -O >> $output1
grep 'nan' $infile | awk '{print $1, $2}' | gmt psxy -R$range -J$projection -O -K -Gdarkblue -Sc0.05 >> $output1


python print_over_value.py $infile 1.5


# THE SECOND MAP
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0weSn -Dh -K -O -X$offset -Y0 >> $output1
gmt pscoast -R$range -J$projection -Lf-124.9/39.15/39.15/50+jt -N1 -N2 -B+t"Average North Offsets" -Wthinner,black -Dh -K -O >> $output1 # the title goes here

# Add data
awk '{print $2, $3, $5}' $infile | gmt psxy -R$range -J$projection -Sc0.3 -Ccolors.cpt -Wthinnest,black -K -O >> $output1
awk '{print $2, $3, $4, $5}' $infile | gmt psvelo -R$range -J$projection -Se0.2/0.68/0 -A+e+gblack+pthickest -Wthick,black -K -O >> $output1
grep 'nan' $infile | awk '{print $1, $2}' | gmt psxy -R$range -J$projection -O -K -Gdarkblue -Sc0.05 >> $output1

gmt psscale -R$range -J$projection -DjTR+w8c/0.5c+o-1.5/1.5 -Ccolors.cpt -B1:"Average Offset":/:mm: -P -O >> $output1
# Convert to PNG
gmt psconvert $output1 -Tg


# Vertical
gmt makecpt -T-4/4/0.2 -Cjet > colors.cpt

# THE FIRST MAP
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output2
gmt pscoast -R$range -J$projection -Lf-124.9/39.15/39.15/50+jt -N1 -N2 -B+t"Average Vert Offsets" -Wthinner,black -Dh -K -O >> $output2 

# Add data
awk '{print $2, $3, $6}' $infile | gmt psxy -R$range -J$projection -Sc0.3 -Ccolors.cpt -Wthinnest,black -K -O >> $output2
awk '{print $2, $3, $1}' $infile | gmt pstext -R$range -J$projection -D0/.15 -F+f6p,Helvetica -K -O >> $output2
awk '{print $2, $3, $4, $5}' $infile | gmt psvelo -R$range -J$projection -Se0.2/0.68/0 -A+e+gblack+pthickest -Wthick,black -K -O >> $output2

gmt psscale -R$range -J$projection -DjTR+w8c/0.5c+o-1.5/1.5 -Ccolors.cpt -B1:"Average Offset":/:mm: -P -O >> $output2




gmt makecpt -T0/2/0.2 -Cjet > colors.cpt

# THE MAP
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0WeSn -Dh -K -X2 -Y2 > $output3
gmt pscoast -R$range -J$projection -Lf-124.9/39.15/39.15/50+jt -N1 -N2 -B+t"Sigma East Offsets" -Wthinner,black -Dh -K -O >> $output3 # the title goes here

# Add data
awk '{print $2, $3, $7}' $infile | gmt psxy -R$range -J$projection -Sc0.3 -Ccolors.cpt -Wthinnest,black -K -O >> $output3
awk '{print $2, $3, $1}' $infile | gmt pstext -R$range -J$projection -D0/.15 -F+f6p,Helvetica -K -O >> $output3
awk '{print $2, $3, $4, $5}' $infile | gmt psvelo -R$range -J$projection -Se0.2/0.68/0 -A+e+gblack+pthickest -Wthick,black -K -O >> $output3
grep 'nan' $infile | awk '{print $1, $2}' | gmt psxy -R$range -J$projection -O -K -Gdarkblue -Sc0.05 >> $output3


# THE SECOND MAP
gmt pscoast -R$range -J$projection -Slightblue -N1 -N2 -B1.0weSn -Dh -K -O -X$offset -Y0 >> $output3
gmt pscoast -R$range -J$projection -Lf-124.9/39.15/39.15/50+jt -N1 -N2 -B+t"Sigma North Offsets" -Wthinner,black -Dh -K -O >> $output3 

# Add data
awk '{print $2, $3, $8}' $infile | gmt psxy -R$range -J$projection -Sc0.3 -Ccolors.cpt -Wthinnest,black -K -O >> $output3
awk '{print $2, $3, $4, $5}' $infile | gmt psvelo -R$range -J$projection -Se0.2/0.68/0 -A+e+gblack+pthickest -Wthick,black -K -O >> $output3
grep 'nan' $infile | awk '{print $1, $2}' | gmt psxy -R$range -J$projection -O -K -Gdarkblue -Sc0.05 >> $output3

gmt psscale -R$range -J$projection -DjTR+w8c/0.5c+o-1.5/1.5 -Ccolors.cpt -B1:"Sigma Offset":/:mm: -P -O >> $output3





open $output3