#!/bin/bash

list=$1

nf=$(head -n 1 $list | awk '{print NF}')

samples=$(ls $list | wc -l)

echo "samples="$samples

echo $list

paste $list > zzz

gawk -v nf="$nf" '{
	for(j=1;j<=nf;j++)
	{
		acum[j]=0.0;
		count[j]=0;
	}; 
	for(i=1;i<=NF;i++)
	{
		acum[(i-1)%nf+1]+=$i;
		count[(i-1)%nf+1]++;
	}; 
	if(NF>0)
	{
		for(j=1;j<=nf;j++){
			printf("%1.2f ",acum[j]/count[j]);
		}
		printf("\n");
	} 
}' zzz

