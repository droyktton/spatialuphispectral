#smo=100000; 
sm(x)=int(x*smo)/(smo+0.); 
set multi lay 2,2; set logs; 


do for [z in "10000000. 1.5 2.0 2.5"] {
	set tit 'z='.z; 
	plot for[i=i0:i1:is] 'Sq_vs_t.dat' index i u (sm($1)*$4**(1./z)):($2/$4**(2.0/z)) w lp smooth un t 'i='.i, 1./x**2
}; 

unset multi
