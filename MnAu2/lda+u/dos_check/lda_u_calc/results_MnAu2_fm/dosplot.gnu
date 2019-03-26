#!/bin/bash

 ef=15.468

 set title "MnAu_2 FM system (LDA+U) - DOS" 
 set xlabel "Energy(Ev)"
 set ylabel "DOS"
 
 set xrange [-15:4]
 
 set arrow 1 nohead from 0,-5 to 0,5 dt 1
 set arrow 2 nohead from -15,0 to 4,0  
 set label 1 "{/Symbol E}f" at 0.015,2

 
 set term png
 set output "mnau2.dos.png"

 plot 'mnau2.dos' u ($1-ef):2 w l notitle, '' u ($1-ef):(-$3) w l notitle
