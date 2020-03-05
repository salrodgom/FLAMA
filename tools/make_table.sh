#!/bin/bash
# 8753.5043 0.418713 1.1082039 : x1 x2 x3
# 7511.9269 0.043775 227.77208 : x4 x5 x6
cc real2binary_string_IEEE.c -lm -o real2binary_string_IEEE
n_refits=0
for x1 in 1.0 0.5 2.0 ; do         
   x01=8800.0
   xx1=$(echo "$x01*$x1" | bc -lq)
   for x2 in 0.90 1.0 1.10 ; do      
   x02=0.4
   xx2=$(echo "$x02*$x2" | bc -lq)
   for x3 in 1.0 0.5 2.0 ; do      
   x03=2.0
   xx3=$(echo "$x03*$x3" | bc -lq)
   #
   for x4 in 1.0 0.5 2.0 ; do      
   x04=7511.9269
   xx4=$(echo "$x04*$x4" | bc -lq)
   #
   for x5 in 0.90 1.0 1.1 ; do      
   x05=0.043775
   xx5=$(echo "$x05*$x5" | bc -lq)
   #
   for x6 in 1.0 0.5 2.0 ; do      
   x06=227.77208
   xx6=$(echo "$x06*$x6" | bc -lq)
   #
   let n_refits++
   xxx1=$(./real2binary_string_IEEE $xx1)
   xxx2=$(./real2binary_string_IEEE $xx2)
   xxx3=$(./real2binary_string_IEEE $xx3)
   xxx4=$(./real2binary_string_IEEE $xx4)
   xxx5=$(./real2binary_string_IEEE $xx5)
   xxx6=$(./real2binary_string_IEEE $xx6)
   echo $xxx1 "#" $xx1 $n_refits
   echo $xxx2 "#" $xx2
   echo $xxx3 "#" $xx3
   echo $xxx4 "#" $xx4
   echo $xxx5 "#" $xx5
   echo $xxx6 "#" $xx6
   done
   done
   done
  done
 done
done
