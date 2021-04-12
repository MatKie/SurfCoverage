#!/bin/bash

begin=5000
cd ..

for dir in $(ls -d VLE_*) 
do
    cd $dir/prod_nvt
    gmx density -f traj.trr -b $begin -dens number -ng 6 -sl 400 -n ../config/index.ndx << EOF
W2
NA+
SO4V9
CM
CT
CM_CT
EOF

cd ../..
done 

for dir in $(ls -d VLE_*) 
do
    cd $dir/prod_nvt 
    gmx energy -f ener.edr -b $begin > energies.out << EOF
1
2
3
4
5
6
7
8
9
10
11
21
25
29
30
32
34
36
38
40
41
42
44
46
47
48
50
52
54
56
58
59
60
EOF

cd ../..

done
