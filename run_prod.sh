#!/bin/bash

for file in $(ls -d VLE_*)
do
    if [ -f $file/prod_nvt/md.log ] 
    then

        if grep -c 'Finished mdrun' $file/prod_nvt/md.log
        then
            echo $file ' done'
            continue
        fi
    fi
    cd $file/prod_nvt
    qsub prod_${file}.sh
    cd ../..
done
