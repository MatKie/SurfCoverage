#!/bin/bash

begin=5000
cd ..
while read line ; 
do
    dir=$line
    cd $dir/prod_vle
    gmx make_ndx -f topol.tpr -o index_analysis.ndx << EOF
del 9
del 8
t W2
t SO4V9
t NA+
t CM | t CT & r SDS
q
EOF

    gmx trjconv -s topol.tpr -f traj.trr -o traj.xtc -pbc mol -center yes  << EOF
SOL
System
EOF
 
    gmx density -f traj.xtc -center yes -b $begin -dens number -ng 4 -sl 400 -n index_analysis.ndx << EOF
W2
W2
NA+
SO4V9
CM_CT_&_SDS
EOF

    rm traj.xtc
    gmx trjconv -f traj.trr -o traj.xtc -skip 3

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

done < vle.done
