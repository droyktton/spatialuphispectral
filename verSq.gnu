unset multi

set logs; 
smo = 5000.
sm(x)=int(x*smo)/(smo+0.); 


set terminal qt 0 font "Sans,9" size 600,700


set multi
set key bot left


#universality1
#is=18
#ie=25

#universality2
is=16
ie=25
istep=1



# sin rescaleo
set yrange [1e-10:1e15]
set xrange [2e-4:pi]
set size 1,0.5
set origin 0,0.5
set title "Sin rescaleo"
set xlabel "q"
set ylabel "S(q,t)"
z=1000000.0; 
plot for[i=1:ie:istep] 'Sq_vs_t.dat' u (t=$4, sm($1)*$4**(1./z)):(t=$4, $2/$4**(2./z)) index i w l t sprintf("t/dt=%d",2**i) smooth un,\
50000./x**2




# con rescaleo z=1.5
z=1.5; 
#set yrange [1e0:1e7]
set yrange [1e0:1e10]
set xrange [1e-2:100]
set size 0.5,0.5
set origin 0,0
set title "Con rescaleo z=1.5"
set xlabel "q t^{1/z}"
set ylabel "S(q,t)/t^{2/z}"
plot for[i=is:ie:istep] 'Sq_vs_t.dat' u (t=$4, sm($1)*$4**(1./z)):(($1<0.01 && $2>1.e5)?($2/$4**(2./z)):(1./0)) \
index i t '' smooth un w lp,\
1000/x**2


# con rescaleo z=2.0
z=2.0;
#set yrange [1e2:1e7]
set yrange [1e2:1e10]
set xrange [1e-2:100]
set size 0.5,0.5
set origin 0.5,0
set title "Con rescaleo z=2.0"
set xlabel "q t^{1/z}"
set ylabel "S(q,t)/t^{2/z}"
plot for[i=is:ie:istep] 'Sq_vs_t.dat' u (t=$4, sm($1)*$4**(1./z)):(($1<0.01 && $2>1.e5)?($2/$4**(2./z)):(1./0)) \
index i t '' smooth un  w lp,\
1000/x**2


unset multi

