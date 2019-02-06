#!/bin/bash


infile1='Fall2018/20140310.txt'
infile2='Fall2018/20161208.txt'
out_ps_name='Figures/combined_2014_2016'
./timing_figure_combined.gmt $infile1 $infile2 -125 -118 36.5 42.0 $out_ps_name # northern California

