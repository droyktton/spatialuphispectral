#!/bin/bash

for((seed=1234;seed<1284;seed++))
do 
	echo $seed 
	./spatialuphispectral 1024 3.0 1048576 $seed 10000000000 0.27
	tar zcvf "L1024_datos"$seed".tar.gz" *.dat
	rm *.dat
done
