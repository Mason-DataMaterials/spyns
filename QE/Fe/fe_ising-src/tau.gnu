#!/bin/bash

set term png
set output 'plots_40/tau.png'

set ylabel '{/Symbol t}'
set xlabel 'T'
set title 'Correlation Time'

 plot 'fit.txt' u 1:3 w lp notitle

