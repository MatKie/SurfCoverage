#!/bin/bash

if [ -f lle.started ];
then 
    rm lle.started
fi
if [ -f lle.done ];
then 
    rm lle.done
fi
if [ -f lle.error ];
then 
    rm lle.error
fi

while read line;
do
    if grep -q 'Statistics over 10001' $line/prod_nvt/md.log
    then
        echo $line >> lle.started
    elif grep -q 'Statistics over 1500001' $line/prod_nvt/md.log
    then 
        echo $line >> lle.done
    else
        echo $line >> lle.error
    fi
done < eqnpt.done
    
