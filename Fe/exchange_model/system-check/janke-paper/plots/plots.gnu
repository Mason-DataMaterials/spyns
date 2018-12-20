#!/bin/bash


set term png

beta(x) = 1.0/x

set grid
#set xrange[1:10]

n12=1728
n16=4096
n20=8000
n24=13824
n32=32768
set output "M.png"
set ylabel "<M>"
set xlabel "T"
plot "12" u 2:3  w l lt 2 t "12",\
     "16" u 2:3  w l  lt 3 t "16",\
     "20" u 2:3  w l  lt 4 t "20",\
     "24" u 2:3  w l  lt 5 t "24",\
     "32" u 2:3  w l  lt 6 t "32",\

     
set output "H.png"
set ylabel "H"
set xlabel "T"
plot "12" u 2:($7/n12)  w l lt 2 t "12",\
     "16" u 2:($7/n16)  w l  lt 3 t "16",\
     "20" u 2:($7/n20)  w l  lt 4 t "20",\
     "24" u 2:($7/n24)  w l  lt 5 t "24",\
     "32" u 2:($7/n32)  w l  lt 6 t "32",\

set output "X.png"
set ylabel "X"
set xlabel "T"
plot "12" u 2:11  w l lt 2 t "12",\
     "16" u 2:11  w l  lt 3 t "16",\
     "20" u 2:11  w l  lt 4 t "20",\
     "24" u 2:11  w l  lt 5 t "24",\
     "32" u 2:11  w l  lt 6 t "32",\

