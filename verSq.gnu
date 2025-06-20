set logs; 
sm(x,f)=int(x*f)/(f+0.); 
z=1000000.0; 

plot for[i=17:27:1] 'Sq_vs_t.dat' u (t=$4, sm($1,1000.)*$4**(1./z)):($2/$4**(2./z)) index i w l t sprintf("t=%1.2f",t) smooth un,\
100/x**2
