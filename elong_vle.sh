#!/bin/bash


while read line;
do
    cd $line
    cd prod_vle
    gmx convert-tpr -s topol.tpr -o topol.tpr -extend 150000
    sed -i s/12/72/g prod_${line}.sh
    cd ../..
done < elong_vle.list
