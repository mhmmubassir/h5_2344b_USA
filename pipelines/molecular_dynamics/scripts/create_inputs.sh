#!/bin/bash
if [ $# -eq 4 ]; then
    if [ -d $1 ]; then
        rm -r $1
    fi
    mkdir $1
    cp $2 $1/structure.pdb
    cp $3 $1/bonding.txt
    cd $1
    cp ../z.template/* .
    cp leapTop.in leap.in
    cat bonding.txt >> leap.in
    cat leapBot.in >> leap.in
    #Leap
    tleap -f leap.in
    #just write out solute
    nla=$((`grep -v TER CPLX.pdb | grep -c "."` -1)) #CPLX.pdb is the file with ATOM and TER card only 
    sed -i 's/ntwprt = 0/ntwprt = '$nla'/g' 10.produ.in #ntpwrt needs the number of atoms to be defined 
    #Jobname
    sed -i 's/ JOBNAME/ '$4'/g' *submit*
    sbatch --no-requeue submit.CPU.sapelo2.sh
    cd ../
else 
    echo "Usage ./create_inputs.sh foldername input.pdb input.link jobname"
fi
