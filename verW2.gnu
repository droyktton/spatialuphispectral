unset multi
set terminal qt 0 font "Sans,9" size 800,900

set key left bottom

set logs

tmax=100000
tmin=5000
wmin=1.0

set multi lay 2,1

set xlabel "t"
set ylabel "w^2(t)"
plot 'z_vs_t.dat' u 1:4 w lp

set xlabel "t"
set ylabel "[w^2(t)-1]/t^{2 {/Symbol z}/z}"
plot [tmin:tmax][0.03:0.4] for[z in "1.5 1.75 2.0"] 'z_vs_t.dat' u 1:(($4>3 && $1<tmax)?(($4-wmin)/$1**(2*0.5/z)):(1/0.)) w lp t 'z='.z

unset multi 
