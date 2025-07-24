#periodico en una islita
./spatialuphispectral 257 3.0 10000000000 1234 10000000000 0.27
#plot  [][:] '< tail -n 100000 bloch.dat' u 1:2:(sin($9)) lc palette ps .5 pt 7


#chaos
./spatialuphispectral 4096 3.0 10000000000 1234 10000000000 0.27
#plot  [][:500] '< tail -n 1000000 bloch.dat' u 1:2:(sin($9)) lc palette ps .1 pt 7

