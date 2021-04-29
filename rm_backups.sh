#!/bin/bash

while read line;
do

    cd $line
    cd eq_npt
    rm \#*
    cd ../prod_nvt
    rm \#*
    cd ../prod_vle
    rm \#*
    cd ../em
    rm \#*
    cd ../config
    rm \#*
    cd ..
    if [ -d eq_npt_hp ];
    then
        cd eq_npt_hp
        rm \#*
        cd ..
    fi
    cd ..

done < list.all
