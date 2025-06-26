#!/bin/bash

list=$1
nf=$2

#nf=$(head -n 1 $list | awk '{print NF}')

samples=$(ls $list | wc -l)

#echo $list 
#echo "nf="$nf" samples="$samples

paste $list > zzz

gawk -v nf="$nf" '{
	for(j=0;j<nf;j++)
	{
		acum[j]=0.0;
		acum2[j]=0.0;
		count[j]=0;
	}; 
	for(i=1;i<=NF;i++)
	{
		acum[(i-1)%nf]+=$i;
		acum2[(i-1)%nf]+=$i*$i;
		count[(i-1)%nf]++;
	}; 
	if(NF>0)
	{
		for(j=1;j<=nf;j++){
			media=acum[j-1]/count[j-1];
			printf("%f ",media);
		}
		for(j=1;j<=nf;j++){
			media=acum[j-1]/count[j-1];
			var=acum2[j-1]/count[j-1]-media*media);
			printf("%f ",var);
		}
		printf("%d \n", count[0]);
	}
	if(NF==0) printf("\n");
}' zzz

