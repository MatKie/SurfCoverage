#!/bin/bash

begin=5000
cd ..
while read line ; 
do
    dir=$line
    cd $dir/prod_nvt
    gmx trjconv -f traj.trr -o traj.xtc -skip 3

    gmx make_ndx -f topol.tpr -o index_analysis.ndx << EOF
del 11
del 10
del 9
t W2
t SO4V9
t NA+
t CM | t CT & r SDS
t CM & r DOD
t CT & r DOD
q
EOF
 
    gmx density -f traj.trr -b $begin -dens number -ng 7 -sl 400 -n index_analysis.ndx << EOF
W2
NA+
SO4V9
CM_CT_&_SDS
DOD
CM_&_DOD
CT_&_DOD
EOF

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

done < lle.done
