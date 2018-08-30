
#!/bin/bash

mcfname='data_10.txt'
mftfname='mft_10.txt'
set term png
set output 'mft.png' 

set xlabel "T"
set key autotitle columnheader
set title "Magnetization"
set ylabel "M"
plot mcfname u 1:($2/(10*10)) w l t 'MC Simulation',     mftfname u 1:2 w p pt 8 t 'MFT Calculation' 



