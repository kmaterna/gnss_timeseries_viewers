#!/bin/bash



range="-125/-121/37.0/42.2"
projection="M3.7i"
out_dila="Strain/delaunay_dilatation.ps"
pdf_out_strain="Strain/delaunay_strain.pdf"
folder="Fields/"


# # Dilatation
gmt makecpt -Iz -D -T-350/350/35 -Cpolar.cpt > mycpt.cpt
# gmt makecpt -Iz -D -T-40/40/2 -Cpolar.cpt > mycpt.cpt
gmt psxy Strain/Dilatation.txt -R$range -J$projection -BWeSN+t"Dilatation" -Bp1.0 -Cmycpt.cpt -Wblack -L -K -P > $out_dila
gmt pscoast -R$range -J$projection -Wthick,black -N2 -Df -Sgray -K -O -P >> $out_dila
gmt psscale -DjTR+w4.5/0.5+o-1.1/1.5 -R$range -J$projection -B200:"Dilatation":/:: -Cmycpt.cpt -K -O -P >> $out_dila
# gmt psvelo positive_eigs.txt -Se0.003/0.68/0 -A+e+pthick,blue -Gblue -R$range -J$projection -K -O >> $out_dila
# gmt psvelo negative_eigs.txt -Se0.003/0.68/0 -A+b+pthick,black -Gred -R$range -J$projection -K -O >> $out_dila


gmt makecpt -T-3/3/0.05 -Cpolar -D > vert.cpt
awk '{print $1, $2, $5}' $folder"MIDAS.txt" | gmt psxy -R$range -J$projection -Sc0.3 -Cvert.cpt -O -K >> $out_dila
gmt psscale -R$range -J$projection -D4.0i/0.4i+w6.5/0.5 -S -B0.5:"Vertical":/:"mm/yr": -Cvert.cpt -K -O >> $out_dila



# Add the fault map (on land)
gmt psxy ../../Reports/Mapping_files/Quaternary.txt -R$range -J$projection -Wthin,gray34 -O -P >> $out_dila
# gmt psvelo tempgps.txt -Se0.03/0.68/0 -A+e+pthick,white -Gwhite -R$range -J$projection -K -O >> $out_strain

gmt psconvert $out_dila -Tf

open $out_dila