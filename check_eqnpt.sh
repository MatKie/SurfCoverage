#!/bin/bash

if [ -f eqnpt.started ];
then 
    rm eqnpt.started
fi
if [ -f eqnpt.done ];
then 
    rm eqnpt.done
fi
if [ -f eqnpt.error ];
then 
    rm eqnpt.error
fi

while read line;
do
    if grep -q 'Statistics over 10001' $line/eq_npt/md.log
    then
        echo $line >> eqnpt.started
    elif grep -q 'Statistics over 2500001' $line/eq_npt/md.log
    then 
        echo $line >> eqnpt.done
    else
        echo $line >> eqnpt.error
    fi
done < started.list
    
