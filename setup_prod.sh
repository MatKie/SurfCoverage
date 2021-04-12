#!/bin/bash

for file in $(ls -d VLE_*)
do
    if [ -f $file/eq_nvt/md.log ] 
    then

        if grep -c 'Finished mdrun' $file/eq_nvt/md.log
        then
            echo $file ' done'
            cp template_VLE/prod_nvt/nvt.mdp $file/prod_nvt
            cd $file/prod_nvt/
            gmx grompp -f nvt.mdp -c ../eq_nvt/confout.gro -p ../topo/topol.top -n ../config/index.ndx
            mv GROMACS.sh prod_${file}.sh
            cd ../..
            continue
        fi
    fi
    echo ${file}' not finished'
done
