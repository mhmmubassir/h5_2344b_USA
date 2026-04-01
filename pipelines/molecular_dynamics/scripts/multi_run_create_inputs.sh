#!/bin/bash

if [ $# -ne 2 ]; then
    echo ""
    echo "Usage: replicates systemList"
    echo "Exmpl: 5 systemList.txt" 
    echo "note names in systemList.txt should correspond to files like y.input_files/\$name.pdb and y.input_files/\$name.bond"
    echo ""
    exit 1
fi

replicates=$1
systemList=$2
pdbFile=y.input_files/

#>folderList.txt

j=0
for name in `cat $systemList`
do
    pdbFile=y.input_files/$name.pdb
    bonding=y.input_files/$name.bond
    for ((i = 0 ; i < $replicates ; i++))
    do
        echo "$j.$i.$name" >> folderList.txt
        ./create_inputs.sh $j.$i.$name $pdbFile $bonding $j.$i.$name
    done
    j=$(($j + 1))
done

