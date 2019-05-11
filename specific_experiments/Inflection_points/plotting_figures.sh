#!/bin/bash


# infile1='Fall2018/20140310.txt'
# infile2='Fall2018/20161208.txt'
type_of_inflection="max_curve"
# type_of_inflection="min_slope"
infile1='Outputs_'$type_of_inflection'/20140310.txt'
infile2='Outputs_'$type_of_inflection'/20161208.txt'
out_ps_name='Figures/combined_2014_2016_'$type_of_inflection
./timing_figure_combined.gmt $infile1 $infile2 -125 -119.4 38.2 42.0 $out_ps_name # northern California

