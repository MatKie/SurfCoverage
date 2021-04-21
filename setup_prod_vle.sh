#!/bin/bash

while read line;
do

    
    if grep -q 'Statistics over' $line/prod_vle/md.log ;
    then
        echo 'Already done' $line
    else
    
    cp -r template_VLE/prod_vle $line
    cd $line
    cd prod_vle
    mv GROMACS.sh prod_$line.sh

    sed s/DOD\ 3400//g ../topo/topol.top > ../topo/topol_vle.top
    
    gmx make_ndx -n ../config/index.ndx -o index.ndx << EOF
3 | 4 | 5
q
EOF

    gmx trjconv -f ../prod_nvt/confout.gro -n index.ndx -o input.gro -s ../prod_nvt/topol.tpr  << EOF
NA_SDS_Water
EOF

    gmx grompp -f dummy.mdp -p ../topo/topol_vle.top -c input.gro 
    gmx make_ndx -f topol.tpr << EOF
t W2
t SO4V9
t NA+
t CM
t CT
t CM | t CT
q
EOF

    gmx grompp -f nvt.mdp -c input.gro -p ../topo/topol_vle.top -n index.ndx
    gmx mdrun -v -nsteps 10000 -s topol.tpr

    cd ../..
    
    fi
done < lle.done
