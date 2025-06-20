./a.out 1.001 1.500 500 0.0 2.0 500

#gawk '{max=($3>$4)?($3):($4); maxx=(max>1.0)?(max):(0.0); if(NF==3){print $2,$1,maxx;}else{print}}' rk4_k_sweep.txt > data.txt; 
#gawk '{max=($3>$4)?($3):($4); maxx=(max>1.0)?(max):(max); if(NF==3){print $1,$2,maxx;}else{print}}' rk4_k_sweep.txt > data.txt; 

gnuplot vermapa.gnu
