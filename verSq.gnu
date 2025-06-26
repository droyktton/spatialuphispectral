unset multi


#file='Sq_vs_t.dat'

set logs; 
smo = 100000.
sm(x)=int(x*smo)/(smo+0.); 
#sm(x)=exp(int(log(x)*smo)/(smo+0.)); 

set terminal qt 0 font "Sans,9" size 700,800


set multi
set key bot left


#universality1
#is=18
#ie=25

#universality2
is=13
ie=22
istep=1



# sin rescaleo
set yrange [1:1e6]
set xrange [2e-4:0.1]
set size 1,0.5
set origin 0,0.5
set title "Sin rescaleo"
set xlabel "q"
set ylabel "S(q,t)"
z=1000000.0; 
#set colorbox
#unset log cb      # ensure color axis (colorbar) is linear

plot for[i=10:ie:istep] file u \
(t=$4, sm($1)*$4**(1./z)):($2/$4**(2./z)) index i \
smooth un w lp lc i t sprintf("t/dt=%d",2**i),\
#for[i=1:ie:istep] file u \
#(t=$4, sm($1)*$4**(1./z)):(t=$4, $2/$4**(2./z)):(sqrt($6/$9)):($4) index i  smooth un t '' w error,\
50000./x**2



qmax=0.007
Sqmin=1e-5
Sqmax=1e5

# con rescaleo z=1.5
z=1.5; 
#set yrange [1e0:1e7]
set yrange [Sqmin:Sqmax]
set xrange [1e-2:100]
set size 0.5,0.5
set origin 0,0
set title "Con rescaleo z=1.5"
set xlabel "q t^{1/z}"
set ylabel "S(q,t)/t^{2/z}"
plot for[i=is:ie:istep] file u (t=$4, sm($1)*$4**(1./z)):(($1<qmax && $2>Sqmin)?($2/$4**(2./z)):(1./0)) \
index i t '' smooth un w lp,\
1000/x**2


# con rescaleo z=2.0
z=2.0;
#set yrange [1e2:1e7]
set yrange [Sqmin:Sqmax]
set xrange [1e-2:100]
set size 0.5,0.5
set origin 0.5,0
set title "Con rescaleo z=2.0"
set xlabel "q t^{1/z}"
set ylabel "S(q,t)/t^{2/z}"
plot for[i=is:ie:istep] file u (t=$4, sm($1)*$4**(1./z)):(($1<qmax && $2>Sqmin)?($2/$4**(2./z)):(1./0)) \
index i t '' smooth un  w lp,\
1000/x**2


unset multi

