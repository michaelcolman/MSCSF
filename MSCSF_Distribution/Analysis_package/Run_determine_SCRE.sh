#!/bin/sh

for i in `seq 1 1 $5`
do
    num=$(printf "%04d" $i)
    echo "${1}${num}"
    ./Determine_SCRE ${1}${num}/CRU.dat $2 $3 $4 
done

