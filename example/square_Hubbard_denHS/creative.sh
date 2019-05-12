#!/bin/bash
for i in "$@" ;do
    mkdir ${i}t
    sed -i "7s/[0-9]\.[0-9][0-9]/$i/" dqmc.in
    make
    mpirun -n 2  ./dqmc.x
#cp *.dat $1
    cp  *.dat "${i}t"
    
done
