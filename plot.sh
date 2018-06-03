#!/bin/bash

gawk -v d=$2 -v w=$3 -v os=$4 'function abs(x) { return (x>=0?x:-x) }
   {
           if($0~/# Contour/) nr=0
           if(nr==int(os+w/2) && (d%2)==0) {a[i]=$1; b[i]=$2; c[i]=$3;}
           if(nr==int(os+w/2)-1 && (d%2)==0) {i++; x = $1; y = $2;}
           if(nr==int(os+w/2)+1 && (d%2)==0) r[i]= 180.0*atan2(y-$2, x-$1)
           if(abs(nr-os-w/2)>w/2 && (d%2)==1) print $0
           nr++
   }
   END {   if(d==0) {
                   for(j=1;j<=i;j++)
                   printf "set label %d \"%g\" at %g, %g centre front rotate by %d tc lt %d \n", j, c[j], a[j], b[j], r[j], j
           }
           if(d==2) {
                   printf "plot \"test.dat\" w ima, \"cont.plt\" index 0 w l lt 1,\\\n"
                   for(j=2;j<i;j++) printf "\"\" index %d w l lt %d,\\\n", j-1, j
                   printf "\"\" index %d w l lt %d\n", i-1, i
           }
   }' $1
