set terminal pngcairo size 600,500
set output 'heatmap.png'

set xlabel 'H'
set ylabel 'k'
set zlabel 'log(floquet multiplier)' 
 
set view map
set pm3d map
set pm3d implicit  # assumes grid structure from x/y but data is unsorted

splot 'rk4_k_sweep.txt' using 2:1:(log(($3>$4)?(($3>1)?($3):(0)):(($4>1)?($4):(0)))) with pm3d notitle

#splot 'data.txt' using 2:1:(log($3)) with pm3d notitle
