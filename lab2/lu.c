#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <math.h>

#define DEFAULT_OMP_THREADS 2

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static int double_equals(double a, double b) {
	if (a == b) {
		return 1;
	}
	double absA = fabs(a);
	double absB = fabs(b);
	double diff = fabs(a - b);
	if (a == 0 || b == 0 || absA + absB < DBL_MIN) {
		return diff < (DBL_EPSILON * DBL_MIN);
	} else {
		return diff / (absA + absB) < DBL_EPSILON;
	}
}

static int lu_omp(double *A, int n) {
	int i, j, k, rows, thr_min, thr_max;
	int thr_id = 0;
	int thr_count;
#pragma omp parallel shared(A,n,thr_count) private(i,j,k,thr_id,rows,thr_min,thr_max)
{
#ifdef _OPENMP
	thr_count = omp_get_num_threads();
#endif

#ifdef _OPENMP
	thr_id = omp_get_thread_num();
#endif
	rows = n / thr_count;
	thr_min = thr_id * rows;
	thr_max = rows == 0 ? -1 : thr_min + rows - 1;

	if (thr_id == thr_count - 1 && (n - (thr_max + 1)) > 0) {
		thr_max = n - 1;
	}
	for (k = 0; k < n - 1; ++k) {
		if (k >= thr_min && k <= thr_max) {
//			for (j = MAX(k + 1, thr_min); j < thr_max; ++j) {
			for (j = k + 1; j < n; ++j) {
				A[k * n + j] /= A[k * n + k];
			}
		}
#pragma omp barrier
		for(i = ((k + 1) > thr_min ? k + 1 : thr_min ); i <= thr_max; ++i) {
			for(j = k + 1; j < n; ++j) {
				A[i * n + j] -= A[i * n + k] * A[k * n + j];
			}
		}
	}
}
	return 0;
}

static inline double det_omp(double *LU, int  n) {
	double d = 1.;
#pragma omp parallel for reduction (*:d) 
	for (int i = 0; i < n; ++i) {
		d *= LU[i * n + i];
	}
	return d;
}

int main(int argc, char **argv) {
	int n;

	if (!scanf("%d", &n)) {
		fprintf(stderr, "Cannot read matrix size");
		return -1;
	}
	int threads_count = (argc == 2) ? atoi(argv[1]) : DEFAULT_OMP_THREADS;
	omp_set_num_threads(threads_count);
	
	double *A = (double *)malloc(n * n * sizeof(double));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			scanf("%lf", A + (n * i + j));
		}
	}
	if (double_equals(A[0], 0.)) {
		printf("Cannot apply LU decomposition\n");
		free(A);
	}

	double tim = omp_get_wtime();
	lu_omp(A, n);
	double det_val = det_omp(A, n);
	tim = omp_get_wtime() - tim;
	printf("detA = %lf\n", det_val);
	printf("%.6f\n", tim);

	free(A);

	return 0;
}
