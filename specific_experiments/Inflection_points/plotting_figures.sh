#!/bin/bash


# infile1='Fall2018/20140310.txt'
# infile2='Fall2018/20161208.txt'
infile1='20140310.txt'
infile2='20161208.txt'
out_ps_name='Figures/combined_2014_2016_curvature'
./timing_figure_combined.gmt $infile1 $infile2 -125 -119.4 38.2 42.0 $out_ps_name # northern California

