#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define I_MAT(R,C,N) ((R) * (N) + (C))
#define SIZE_MAT(N) ((N) * (N) * sizeof(double))
#define SIZE_ROW(N) ((N) * sizeof(double))

int my_log2(int x) {
	int result = 1;
	int tmp_x = x;
	while (tmp_x >>= 1) {
		result++;
	}
	if (x == 1 << (result - 1)) {
		return result - 1;
	}
	return result;
}

int get_new_dimension(int x) {
    return 1 << my_log2(x);
}

void multiply(int n, double *a, double *b, double *c) {
	double col_b[n];

	for (int j = 0; j < n; ++j) {
		for (int k = 0; k < n; ++k) {
			col_b[k] = b[I_MAT(k,j,n)];
		}
		for (int i = 0; i < n; ++i) {
			double *row_a = a + i * n;
			double sum = 0.;
			for (int k = 0; k < n; ++k) {
				sum += row_a[k] * col_b[k];
			}
			c[I_MAT(i,j,n)] = sum;
		}
	}

	/*for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double sum = 0;
			for (int k = 0; k < n; ++k) {
				sum += a[I_MAT(i,k,n)] * b[I_MAT(k,j,n)];
			}
			c[I_MAT(i,j,n)] = sum;
		}
	}*/
}

void summation(int n, double *a, double *b, double *c) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			c[I_MAT(i,j,n)] = a[I_MAT(i,j,n)] + b[I_MAT(i,j,n)];
		}
	}
}

void subtraction(int n, double *a, double *b, double *c) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			c[I_MAT(i,j,n)] = a[I_MAT(i,j,n)] - b[I_MAT(i,j,n)];
		}
	}
}

void split_matrix(int n, double *m, double *m11, double *m12, double *m21, double *m22) {
	int new_n = n >> 1;
	for (int i = 0; i < new_n; ++i) {
		memcpy(m11 + i * new_n, m + i * n, SIZE_ROW(new_n));
		memcpy(m12 + i * new_n, m + i * n + new_n, SIZE_ROW(new_n));
		memcpy(m21 + i * new_n, m + (i + new_n) * n, SIZE_ROW(new_n));
		memcpy(m22 + i * new_n, m + (i + new_n) * n + new_n, SIZE_ROW(new_n));
	}
}

void collect_matrix(int n, double *m, double *m11, double *m12, double *m21, double *m22) {
	int new_n = n << 1;
	for (int i = 0; i < n; ++i) {
		memcpy(m + i * new_n, m11 + i * n, SIZE_ROW(n));
		memcpy(m + i * new_n + n, m12 + i * n, SIZE_ROW(n));
		memcpy(m + (i + n) * new_n, m21 + i * n, SIZE_ROW(n));
		memcpy(m + (i + n) * new_n + n, m22 + i * n, SIZE_ROW(n));
	}
}

int  multi_strassen(int n, double *a, double *b, double *c) {
	if (n <=128) {
		multiply(n, a, b, c);
		return 0;	
	}
	n >>= 1;

	double *matrices = (double *)malloc((4 + 5 + 4 + 5 + 7 + 4 + 4) * SIZE_MAT(n));
	if (!matrices) {
		return 1;
	}
	double *a11 = matrices;
	double *a12 = a11 + n * n;
	double *a21 = a12 + n * n;
	double *a22 = a21 + n * n;

	double *a11_sum_a12 = a22 + n * n;
	double *a11_sum_a22 = a11_sum_a12 + n * n;
	double *a21_sum_a22 = a11_sum_a22 + n * n;
	double *a12_sub_a22 = a21_sum_a22 + n * n;
	double *a21_sub_a11 = a12_sub_a22 + n * n;
	
	double *b11 = a21_sub_a11 + n * n;
	double *b12 = b11 + n * n;
	double *b21 = b12 + n * n;
	double *b22 = b21 + n * n;

	double *b11_sum_b12 = b22 + n * n;
	double *b11_sum_b22 = b11_sum_b12 + n * n;
	double *b21_sum_b22 = b11_sum_b22 + n * n;
	double *b12_sub_b22 = b21_sum_b22 + n * n;
	double *b21_sub_b11 = b12_sub_b22 + n * n;

	double *p1 = b21_sub_b11 + n * n;
	double *p2 = p1 + n * n;
	double *p3 = p2 + n * n;
	double *p4 = p3 + n * n;
	double *p5 = p4 + n * n;
	double *p6 = p5 + n * n;
	double *p7 = p6 + n * n;

	double *p1_sum_p4 = p7 + n * n;
	double *p7_sub_p5 = p1_sum_p4 + n * n;
	double *p1_sub_p2 = p7_sub_p5 + n * n;
	double *p3_sum_p6 = p1_sub_p2 + n * n;
	double *c11 = p3_sum_p6 + n * n;
	double *c12 = c11 + n * n;
	double *c21 = c12 + n * n;
	double *c22 = c21 + n * n;

	split_matrix(n * 2, a, a11, a12, a21, a22);
	split_matrix(n * 2, b, b11, b12, b21, b22);

	summation(n, a11, a22, a11_sum_a22);
	summation(n, b11, b22, b11_sum_b22);
	summation(n, a21, a22, a21_sum_a22);
	subtraction(n, b12, b22, b12_sub_b22);
	subtraction(n, b21, b11, b21_sub_b11);
	summation(n, a11, a12, a11_sum_a12);
	subtraction(n, a21, a11, a21_sub_a11);
	summation(n, b11, b12, b11_sum_b12);
	subtraction(n, a12, a22, a12_sub_a22);
	summation(n, b21, b22, b21_sum_b22);

	if (multi_strassen(n, a11_sum_a22, b11_sum_b22, p1)
			|| multi_strassen(n, a21_sum_a22, b11, p2)
			|| multi_strassen(n, a11, b12_sub_b22, p3)
			|| multi_strassen(n, a22, b21_sub_b11, p4)
			|| multi_strassen(n, a11_sum_a12, b22, p5)
			|| multi_strassen(n, a21_sub_a11, b11_sum_b12, p6)
			|| multi_strassen(n, a12_sub_a22, b21_sum_b22, p7)) {
		free(matrices);
		return 1;
	}

	summation(n, p1, p4, p1_sum_p4);
	subtraction(n, p7, p5, p7_sub_p5);
	subtraction(n, p1, p2, p1_sub_p2);
	summation(n, p3, p6, p3_sum_p6);

	summation(n, p1_sum_p4, p7_sub_p5, c11);
	summation(n, p3, p5, c12);
	summation(n, p2, p4, c21);
	summation(n, p1_sub_p2, p3_sum_p6, c22);

	collect_matrix(n, c, c11, c12, c21, c22);

	free(matrices);
	return 0;
}

int main(int argc, char **argv) {
	int m_dim_nopad = atoi(argv[argc - 1]);
	int m_dim = get_new_dimension(m_dim_nopad);
	double *a = (double *)malloc(SIZE_MAT(m_dim));
	if (!a) {
		goto errA;
	}
	double *b = (double *)malloc(SIZE_MAT(m_dim));
	if (!b) {
		goto errB;
	}
	double *c = (double *)malloc(SIZE_MAT(m_dim));
	if (!c) {
		goto errC;
	}


	/*for (int i = 0; i < m_dim; ++i) {
		for (int j = 0; j < m_dim; ++j) {
			if (i < m_dim_nopad && j < m_dim_nopad) {
				scanf("%lf", &(a[I_MAT(i,j,m_dim)]));
			} else {
				a[I_MAT(i,j,m_dim)] = 0.;
			}
		}
	}
	for (int i = 0; i < m_dim; ++i) {
		for (int j = 0; j < m_dim; ++j) {
			if (i < m_dim_nopad && j < m_dim_nopad) {
				scanf("%lf", &(b[I_MAT(i,j,m_dim)]));
			} else {
				b[I_MAT(i,j,m_dim)] = 0.;
			}
			c[I_MAT(i,j,m_dim)] = 0.;
		}
	}*/
	for (int i = 0; i < m_dim; ++i) {
		for (int j = 0; j < m_dim; ++j) {
			if (i == j && i < m_dim_nopad) {
				a[I_MAT(i,j,m_dim)] = 1.;
				b[I_MAT(i,j,m_dim)] = 1.;
			} else {
				a[I_MAT(i,j,m_dim)] = 0.;
				b[I_MAT(i,j,m_dim)] = 0.;
			}
			c[i * m_dim + j] = 0.;
		}
	}

	clock_t t = clock();
	if (multi_strassen(m_dim, a, b, c)) {
		printf("Cannot eval C=A*B\n");
		goto errC;
	}
	t = clock() - t;
	double time_taken = ((double)t / CLOCKS_PER_SEC);
	printf("Execution time: %.6f\n", time_taken);


	/*for (int i = 0; i < m_dim_nopad; ++i) {
		for (int j = 0; j < m_dim_nopad; ++j) {
			printf("%lf ", c[i * m_dim + j]);
		}
		printf("\n");
	}*/

	free(a);
	free(b);
	free(c);

	return 0;
errC: {
	      free(c);
      }
errB: {
	      free(b);
      }
errA: {
	      free(a);
	      return 1;
      }
}
