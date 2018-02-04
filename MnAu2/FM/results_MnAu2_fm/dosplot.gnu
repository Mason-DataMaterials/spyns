#!/usr/bin/gnuplot -persist
set terminal eps
set output 'mnau2_fm-dos.eps'

set label 1 "{/Symbol E}_f" at 0.00000, 5.40000, 0.00000 left norotate back nopoint

set arrow 1 from 0.00000, -6.00000, 0.00000 to 0.00000, 6.00000, 0.00000 nohead back nofilled linecolor rgb "dark-violet"  linewidth 1.500 dashtype 2
set tics scale  1, 0.5, 1, 1, 1
set title "MnAu_2 DOS -- NM" 
set xlabel "Energy (Ev)" 
set xrange [ -20.0000 : 20.0000 ] noreverse nowriteback
set ylabel "DOS" 
set yrange [ -6.00000 : 6.00000 ] noreverse nowriteback
EF = 15.317
x = 0.0
plot "mnau2.pdos_tot" u ($1-EF):2 w l t 'Up', "" u ($1-EF):(-$3) w l t 'Down'

#    EOF
