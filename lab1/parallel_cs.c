#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv) {
	int n = 1000000000;
	double h = 1. / (double)n;
	double sum = 0.;
	omp_set_num_threads(atoi(argv[1]));
	double tim = omp_get_wtime();
#pragma omp parallel 
{
	int numt = omp_get_num_threads();
#pragma omp for 
	for (int i = 0; i < n; ++i) {
		double x = h * ((double)i - 0.5);
#pragma omp critical
		sum += (4. / (1. + x * x));
	}
}
	double pi = h * sum;
	tim = omp_get_wtime() - tim;
	printf("%.16f,%.6f\n", pi, tim);
	return 0;
}

