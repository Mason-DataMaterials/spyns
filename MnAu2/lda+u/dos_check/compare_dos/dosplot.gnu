#!/bin/bash

 ef=15.319
 ef_lda=15.468

 set title "Comparing MnAu_2 DOS with and without LDA+U" 
 set xlabel "Energy(Ev)"
 set ylabel "DOS(E)"
 
 #set key left 


 set xzeroaxis
 set xrange [-10:5]

 
 set arrow 1 nohead from 0,-5 to 0,5  linecolor rgb "#009e73"  linewidth 1.500 dashtype 4 
 set label 1 "{/Symbol E}f" at 0.015,2

 
 set term png
 set output "mnau2.dos.png"

 plot 'mnau2.dos' u ($1-ef):2 w l lt 2 t "dos-up", '' u ($1-ef):(-$3) w l  lt 2 t "dos-down",\
      'mnau2_lda_u.dos' u ($1-ef_lda):2 w l lt 4 t "dos(lda+u)-up", '' u ($1-ef_lda):(-$3) w l  lt 4 t "dos(lda+u)-down"


 set term eps

 set output "mnau2.dos.eps"
 replot

