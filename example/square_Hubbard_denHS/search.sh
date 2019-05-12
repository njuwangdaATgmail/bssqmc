#!/bin/bash
echo  "spin:"
for i in "$@" ;do 
    echo "$i"
    grep "0.000     0     0     0" $i/spin_r.dat
done
echo  "doublon:"
for i in "$@" ;do 
    echo "$i"
    grep "0.000     0     0     0" $i/doublon_r.dat
done

