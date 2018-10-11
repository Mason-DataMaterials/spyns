#!/bin/bash

#M
fname='data_500/data_500.txt'
set term png
set output 'Results/ising_M.png' 

set xlabel "T"
set key autotitle columnheader
set title "Magnetization(N=500)"
set ylabel "<M>"
set xrange [0:20]
plot fname u 1:2 w lp pt 7 notitle, "" u 1:2:(sqrt($3)) w errorbars notitle



#E
reset
fname='data_500/data_500.txt'
set term png
set output 'Results/ising_E.png' 

set xrange [0:20]
set xlabel "T"
set key autotitle columnheader
set title "Energy(N=500)"
set ylabel "<E>" 
plot fname u 1:6 w lp pt 7 notitle, "" u 1:6:(sqrt($7)) w errorbars notitle

#X
reset
fname='data_500/data_500.txt'
set term png
set output 'Results/ising_X.png' 

set xrange [0:20]
set xlabel "T"
set key autotitle columnheader
set title "Susceptibility(N=500)"
set ylabel "{/Symbol C}" 
plot fname u 1:9 w lp pt 7 notitle

#Cv
reset
fname='data_500/data_500.txt'
set term png
set output 'Results/ising_Cv.png' 

set xlabel "T"

set xrange [0:20]
set key autotitle columnheader
set title "Heat Capacity(N=500)"

set ylabel "C_v" 
plot fname u 1:($10) w lp pt 7 notitle







