paste Sq_vs_t*.dat | \
awk '{
	if(NF>4){
		Su=Sphi=n=0;
		for(i=2;i<=NF;i+=4){
			Su+=$i;Sphi+=$(i+1);n++;
		}; 
		print $1, Su/n, Sphi/n, $4}; 
	if(NF==0) printf("\n\n");
}' > zzz


gnuplot -p -e '
	set term wxt size 1000,400; 
	zeta=1./2; 
	set multi lay 1,3; 
	do for[z in "1.5 2.0 10000000."]{
		plot [][:] "Sq_vs_t.dat" u (q=int(log($1*$4**(1./z))*5)/5.,q):(Sq=log($2/$4**(1./z+2.*zeta/z)), Sq) index 13:20 w lp smooth un, -2*x+8.5;
	}; 
	unset multi'
