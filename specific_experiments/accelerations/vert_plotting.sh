#!/bin/bash

outdir="gldas_lssq_ITRF/"
name="2016"
data_file=$outdir$name.txt

./accel_map_vert.gmt $data_file -124.6 -118.4 35.5 43.2 $outdir"NC"$name