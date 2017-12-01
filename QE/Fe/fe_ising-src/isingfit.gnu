#!/bin/bash

reset
set term png
 T=20
 filename='aucf20.txt'
 outfile='fit20.png'
set output outfile

set xrange [0:20]

f(x) = A0*exp(-x/tau)
A0=1000;tau=1

fit f(x) filename  using 1:($2/$3)  via A0, tau
plot filename u 1:3 w p t '@T=$T', f(x) w l t 'Fit'

set print "fit.txt" append
print  T ,"\t", A0, "\t", tau
