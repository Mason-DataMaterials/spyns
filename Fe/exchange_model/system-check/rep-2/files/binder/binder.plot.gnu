#!/bin/bash
set term png
set output "binder.png"
set xlabel "T(K)"
set ylabel "U(T)"

binder(a,b) = 1.0 -  ( a /3.0*b*b )

set xrange[350:1150]
plot 'm.1024' u 1:( binder($6,$4) ) w l dt 2 t "N=1024",     'm.1458' u 1:( binder($6,$4) ) w l dt 3 t "N=1458",     'm.2000' u 1:( binder($6,$4) ) w l dt 4 t "N=2000",     'm.3456' u 1:( binder($6,$4) ) w l dt 5 t "N=3456"

