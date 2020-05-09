#!/bin/bash
echo "physically_constrained
n_par 18
xxxxxxxxx1  A_buck
xxxxxxxxx2  rho_buck
xxxxxxxxx3  C_buck
xxxxxxxxx4  A_buck
xxxxxxxxx5  rho_buck
xxxxxxxxx6  C_buck
xxxxxxxxx7  A_buck
xxxxxxxxx8  rho_buck
xxxxxxxxx9  C_buck
xxxxxxxx10  A_buck
xxxxxxxx11  rho_buck
xxxxxxxx12  C_buck
xxxxxxxx13  A_buck
xxxxxxxx14  rho_buck
xxxxxxxx15  C_buck
xxxxxxxx16  A_buck
xxxxxxxx17  rho_buck
xxxxxxxx18  C_buck
ffit?  .true.
refit? .true.
1" > input 
# 
echo "8753.5043 
0.418713 
1.1082039 
7511.9269 
0.043775 
227.77208 
8955.2901 
0.303895 
267.49087 
62088.839 
0.275171 
1022.4065 
40322.272 
0.330972 
0.0       
0.0       
0.250000 
1024.359" > parameters 
cat parameters | while read line ; do y=$(./real2binary_string_IEEE $line) ; echo "$y # $line" ; done >> input
rm -rf parameters
