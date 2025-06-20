awk -v n=100 'BEGIN{i=0;s=0;}{a[i%n]+=$1;i++;if(s%n==0) s++;}END{for(i=0;i<n;i++) print a[i]/s;}' zzz.dat > zzz2.dat
