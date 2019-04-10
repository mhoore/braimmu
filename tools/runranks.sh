#!/bin/bash

if [ -z "$EXE" ]; then
	if [ -f "$1" ]; then
		EXE="$1"
	else
		EXE="/p/home/jusers/kelling1/juron/checkout/braimmu/src_mpi/braimmu.exe"
	fi
fi

CORES_PER_SOCKET=10

while read LINE; do
	set $LINE
	sed -e "s/@RX@/$1/" -e "s/@RY@/$2/" -e "s/@RZ@/$3/" < in.tpl > in.run
	LOG=run_${1}_${2}_${3}_${4}.log
	echo "running $LOG"
	NPER=""
	if [ "$4" -le "$CORES_PER_SOCKET" ]; then
		NPER="--npersocket $CORES_PER_SOCKET"
	fi
	echo $( which mpirun ) $NPER -np $4 $EXE in.run | sh 2>&1 | tee $LOG | grep -e "^Step" -e "^Set" -e "^Read"
	echo "exitcode: $?"
done < ranks.dat
