#!/bin/bash

#PBS -l nodes=4:gold6128:ram192gb:ppn=2
cd "$PBS_O_WORKDIR" || exit
echo Launching the parallel job from mother superior "$(hostname)"...
make clean
make
for i in {1..5}
do
	mpirun -machinefile "$PBS_NODEFILE" -n 4 ./build/cannon 2048
done
