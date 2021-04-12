#!/bin/bash

for file in $(ls -d VLE_*)
do
    if [ -f $file/eq_nvt/md.log ] 
    then

        if grep -c 'Finished mdrun' $file/eq_nvt/md.log
        then
            echo $file ' done'
            continue
        fi
    fi
    cd $file/eq_nvt
    qsub SC_${file}.sh
    cd ../..
done
