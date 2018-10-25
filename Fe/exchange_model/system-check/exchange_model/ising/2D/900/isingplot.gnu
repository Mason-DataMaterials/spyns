#!/bin/bash

fname="data_900.txt" 
#M
set term png
set output 'Results/ising_M.png' 

set xlabel "T"
set key autotitle columnheader
set title "Magnetization(N=900)"
set ylabel "<M>"
set xrange [0:15]
plot fname u 1:2 w lp pt 7 notitle, "" u 1:2:(sqrt($3)) w errorbars notitle



#E
set term png
set output 'Results/ising_E.png' 

set xrange [0:15]
set xlabel "T"
set key autotitle columnheader
set title "Energy(N=900)"
set ylabel "<E>" 
plot fname u 1:6 w lp pt 7 notitle, "" u 1:6:(sqrt($7)) w errorbars notitle

#X
set term png
set output 'Results/ising_X.png' 

set xrange [0:15]
set xlabel "T"
set key autotitle columnheader
set title "Susceptibility(N=900)"
set ylabel "{/Symbol C}" 
plot fname u 1:9 w lp pt 7 notitle

#Cv
set term png
set output 'Results/ising_Cv.png' 

set xlabel "T"

set xrange [0:15]
set key autotitle columnheader
set title "Heat Capacity(N=900)"

set ylabel "C_v" 
plot fname u 1:10 w lp pt 7 notitle







