#!/bin/bash

avg() {
	n=1
	awk "function isnum(x){return(x==x+0)} { if(isnum(\$$n)) { sum+=\$$n; sumsq+=\$$n*\$$n ; n+=1;} } END { print sum/n, sum, n, sqrt(sumsq/n - sum*sum/n/n) }"
}

for i in $( find -name 'run_*.log' ); do
	echo $i | sed -e 's/.*run_//' -e 's/\.log//' | tr '_' '\t' | tr '\n' '\t'
	grep '^Step' $i | grep -v 'Step 10' | cut -d ' ' -f 10 | awk '{print 1/$1;}' | avg
done | sort -k 4 -g
