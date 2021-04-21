while read line;
do
    split=$(echo $line | cut -d "_" -f 2)
    ./start_VLE.sh $split
done < eqnpt.error
