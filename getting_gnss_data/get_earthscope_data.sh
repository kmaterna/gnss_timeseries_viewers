#!/bin/bash



# UPDATE CWU/PBO POS, OFFSET, EVENT, AND VELOCITY DATA
# First, enable your es sso token with 'es sso login' (expires every day)
# THIS WILL TAKE A WHILE, and it will look like it's not working at first, but it will be working.
# It is a breadth first search
mkdir -p Earthscope_Data/
mkdir -p Earthscope_Data/Time_Series/
cd Earthscope_Data/Time_Series/
wget --recursive --no-directories --no-parent -N --accept="*.cwu.final_nam14.pos,*.cwu.final_igs14.pos" -e robots=off --level=2 --verbose --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/position/
cd ../../

mkdir -p Earthscope_Data/Earthscope_Event_Files/
cd Earthscope_Data/Earthscope_Event_Files/
wget -r -np --reject=tmp,ps,index* --accept "*_kalts.evt" -e robots=off --no-directories --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/event/
cp data.unavco.org/archive/gnss/products/event/*.evt .
rm -r data.unavco.org
cd ../../

mkdir -p Earthscope_Data/Offsets/
cd Earthscope_Data/Offsets/
wget -N --recursive --no-parent --no-directories --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/offset/cwu.kalts_nam14.off
cd ../../

mkdir -p Earthscope_Data/Velocities/
cd Earthscope_Data/Velocities/
wget -N --recursive --no-parent -e robots=off --no-directories --accept "cwu.final_nam14.vel, cwu.final_igs14.vel" --header "authorization: Bearer $(es sso access --token)" https://data.unavco.org/archive/gnss/products/velocity/
cd ../../
