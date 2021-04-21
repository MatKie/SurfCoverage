#!/bin/bash

while read line;
do
    if grep -q 'Statistics over' $line/prod_nvt/md.log ;
    then
        echo 'Already done' $line
    else

    cp template_VLE/prod_nvt/* $line/prod_nvt 
    cd $line/eq_npt
    gmx energy -f ener.edr -b 5000 > energies.out<< EOF
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
12
13
14
15
16
17
18
28
32
36
37
EOF
    gmx energy -f ener.edr -o full_energy.xvg > full_energies.out<< EOF
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
12
13
14
15
16
17
18
28
32
36
37
39
41
43
45
47
48
49
51
53
54
55
57
59
61
63
65
66
67
EOF

cd ../prod_nvt
python change_boxsize.py
mv GROMACS.sh prod_nvt_$file.sh
gmx grompp -f nvt.mdp -c input.gro -n ../config/index.ndx -p ../topo/topol.top
gmx mdrun -v -nsteps 10000 -s topol.tpr
cd ../..
fi
done < eqnpt.done
