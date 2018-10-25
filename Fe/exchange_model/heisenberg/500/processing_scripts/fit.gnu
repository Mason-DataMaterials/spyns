#!/usr/bin/gnuplot -persist
# set terminal qt 0 font "Sans,9"
# set output
set tics scale  1, 0.5, 1, 1, 1
set title "" 
set xlabel "Time Steps" 
#set xrange [ 0.00000 : 30.0000 ] noreverse nowriteback
set ylabel "{/Symbol C}(t)" 

f(x) = A*exp(-x/tau)
A = 1
tau = 1
T=10 
fit f(x) 'aucf10.txt' u 1:3 via A, tau 

set term png
set output 'fit10.png' 

plot 'aucf10.txt' u 1:3 w p pt 7 t 'T=10', f(x) w l t 'Fit'

set print "fit.txt" append
print  T ,"\t", A, "\t", tau

#    EOF