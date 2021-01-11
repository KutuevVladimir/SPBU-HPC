#!/usr/bin/env bash

gcc -O3 -ffast-math consistent.c -o consistent
gcc -fopenmp parallel.c -o parallel
gcc -fopenmp -O3 parallel.c -o parallel_O3
gcc -fopenmp -O3 -ffast-math parallel.c -o parallel_O3_fm
gcc -fopenmp parallel_cs.c -o parallel_cs
gcc -fopenmp -O3 parallel_cs.c -o parallel_cs_O3
gcc -fopenmp -O3 -ffast-math parallel_cs.c -o parallel_cs_O3_fm

for th_count in {1..8}
do
	echo "Threads count $th_count"
	for i in {1..2} 
	do
		./parallel "${th_count}" >> "parallel_${th_count}.txt"
		./parallel_O3 "${th_count}" >> "parallel_O3_${th_count}.txt"
		./parallel_O3_fm "${th_count}" >> "parallel_O3_fm_${th_count}.txt"
		./parallel_cs "${th_count}" >> "parallel_cs_${th_count}.txt"
		./parallel_cs_O3 "${th_count}" >> "parallel_cs_O3_${th_count}.txt"
		./parallel_cs_O3_fm "${th_count}" >> "parallel_cs_O3_fm_${th_count}.txt"
	done
done

