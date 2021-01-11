#!/bin/bash

#PBS -l nodes=1:gold6128:ram192gb:ppn=2
cd "$PBS_O_WORKDIR" || exit
echo Job from mother superior "$(hostname)"...
make clean
make
for i in {1..5}
do
	./build/strassen 2048
done

