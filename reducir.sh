if [ -f "bloch_lines.dat" ]; then
    rm bloch_lines.dat
fi
if [ -f "numberofbl.dat" ]; then
    rm numberofbl.dat
fi


for((t=0;t<=1000000;t+=1000))
do


	gawk -v t=$t 'BEGIN{prev=0;pi=3.14159;}{if($1==t){piso = int($3/pi+0.5)*pi; if(prev!=piso) print t,$2,$3,$4; prev=piso;}}' bloch.dat > zzz 

	echo $t $(wc -l zzz | awk '{print $1}')
	echo $t $(wc -l zzz | awk '{print $1}') >> numberofbl.dat

	if [ -s zzz ]; then
		cat zzz >> bloch_lines.dat
		printf "\n" >> bloch_lines.dat
	fi
done


