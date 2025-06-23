#!/bin/bash

list=$1

nf=$(head -n 1 $list | awk '{print NF}')

samples=$(ls $list | wc -l)

echo "samples="$samples

paste $list | awk -v nf=$nf '{\
	for(j=0;j<nf;j++)
	{
		acum[j]=0;
	}; 
	for(i=0;i<NF;i++){
		acum[i%nf]+=$i
	}; 
	if(NF>0){
		print acum*1.0/NF; else print;}' 
	} 
> "sofq_"$samples"samples.dat"

file="sofq_"$samples"samples.dat"
echo $file

#gnuplot -p -e "set term png; set out 'sofq.png';set logs; plot for[i=0:7] '$file' index i u 0:1 w lp t ''"
