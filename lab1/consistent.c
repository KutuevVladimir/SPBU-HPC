#include <stdio.h>
#include <time.h>

int main() {
	int n = 1000000000;
	double h = 1. / (double)n;
	double sum = 0.;
	clock_t start_time = clock();
	for (int i = 0; i < n; ++i) {
		double x = h * ((double)i - 0.5);
		sum += (4. / (1. + x * x));
	}
	double pi = h * sum;
	double run_time =((double)(clock() -start_time)) / CLOCKS_PER_SEC;
	printf("%.16f,%.6f\n", pi, run_time);
	return 0;
}

