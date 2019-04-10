#!/bin/bash
for i in {1..7}; do
	for j in $( seq 1 $i ); do
		echo $j $i | awk '{n=$1*$2; if(n < 40) {print $1,$2,1,n;}}'
	done
done | sort -t ' ' -k 4 -g > ranks.dat
