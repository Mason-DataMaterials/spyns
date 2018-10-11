#!/bin/bash

#M
fname='data.txt'
set term png
set output 'Results/heisenberg_r_M.png' 

set xlabel "T"
set key autotitle columnheader
set title "Magnetization(N=1024)"
set ylabel "<M>"
set xrange [1000:1350]
plot fname u 1:2 w lp pt 7 notitle, "" u 1:2:(sqrt($3)) w errorbars notitle



#E
reset
fname='data.txt'
set term png
set output 'Results/heisenberg_r_E.png' 

set xrange [1000:1350]
set xlabel "T"
set key autotitle columnheader
set title "Energy(N=1024)"
set ylabel "<E>" 
plot fname u 1:6 w lp pt 7 notitle, "" u 1:6:(sqrt($7)) w errorbars notitle

#X
reset
fname='data.txt'
set term png
set output 'Results/heisenberg_r_X.png' 

set xrange [1000:1350]
set xlabel "T"
set key autotitle columnheader
set title "Susceptibility(N=1024)"
set ylabel "{/Symbol C}" 
plot fname u 1:9 w lp pt 7 notitle

#Cv
reset
fname='data.txt'
set term png
set output 'Results/heisenberg_r_Cv.png' 

set xlabel "T"

set xrange [1000:1350]
set key autotitle columnheader
set title "Heat Capacity(N=1024)"

set ylabel "C_v" 
plot fname u 1:10 w lp pt 7 notitle







