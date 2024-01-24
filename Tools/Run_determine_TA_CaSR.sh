#!/bin/bash
for l in v2_Orig_Control v2_Orig_Fib_10_25 v2_Orig_Fib_10_50 v2_Orig_Fib_50_25 v2_Orig_Fib_50_50 v2_Orig_Fib_10_25_IK1 v2_Orig_Fib_10_50_IK1 v2_Orig_Fib_50_25_IK1 v2_Orig_Fib_50_50_IK1
do
	for i in `seq 1050 5 1300`
    do
    	num1=$(printf "%04d" $i) # for CaSR
		./Determine_TA_CaSR $l 2000 500 10 $num1
	done
done
