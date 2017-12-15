#!/bin/bash

#M
fname='data_8192.txt'
set term png
set output './Results/plots_8192/heisenberg_M.png' 

set xlabel "T"
set key autotitle columnheader
set title "Magnetization"
set ylabel "<M>"
plot fname u 1:2 w lp pt 7 notitle, "" u 1:2:(sqrt($3)) w errorbars notitle



#E
reset
fname='data_8192.txt'
set term png
set output './Results/plots_8192/heisenberg_E.png' 

set xlabel "T"
set key autotitle columnheader
set title "Energy"
set ylabel "<E>" 
plot fname u 1:6 w lp pt 7 notitle, "" u 1:6:(sqrt($7)) w errorbars notitle

#X
reset
fname='data_8192.txt'
set term png
set output './Results/plots_8192/heisenberg_X.png' 

set xlabel "T"
set key autotitle columnheader
set title "Susceptibility"
set ylabel "{/Symbol C}" 
plot fname u 1:9 w lp pt 7 notitle

#Cv
reset
fname='data_8192.txt'
set term png
set output './Results/plots_8192/heisenberg_Cv.png' 

set xlabel "T"

set key autotitle columnheader
set title "Heat Capacity"

set ylabel "C_v" 
plot fname u 1:10 w lp pt 7 notitle







