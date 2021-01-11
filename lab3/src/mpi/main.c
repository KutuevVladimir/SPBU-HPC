#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

#define I_PROC(R,C) ((R) * p_dim + (C))
#define I_GMAT(R,C) ((R) * m_dim + (C))
#define I_SMAT(R,C) ((R) * sub_dim + (C))

int main(int argc, char **argv) {
	int rank, size;
	double start, end;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int p_dim = sqrt(size);
	if (p_dim * p_dim != size) {
		printf("The number of process must be a square.\n");
		goto main_err;
	}
	int my_row, my_col;
	int *proc = (int *)malloc(p_dim * p_dim * sizeof(int));
	int p = 0;
	for (int i = 0; i < p_dim; ++i){
		for (int j = 0; j < p_dim; ++j){
			proc[I_PROC(i,j)] = p;
			if (p == rank){ // make sure every process knows where in the matrix it is
				my_row = i;
				my_col = j;
			}
			p++;
		}
	}

	int m_dim_nopad = atoi(argv[argc - 1]);
	int m_dim =  p_dim * (m_dim_nopad / p_dim + (m_dim_nopad % p_dim != 0));
	double *full_A;
	double *full_B;
	double *full_C;

	MPI_Alloc_mem(m_dim*m_dim * sizeof(double), MPI_INFO_NULL, &full_A);
	MPI_Alloc_mem(m_dim*m_dim * sizeof(double), MPI_INFO_NULL, &full_B);
	MPI_Alloc_mem(m_dim*m_dim * sizeof(double), MPI_INFO_NULL, &full_C);

	// divide matrices into subblocks and subblocks over processes
	int sub_dim = m_dim / p_dim;
	double *loc_A;
	double *loc_B;
	double *loc_C;
	size_t loc_mat_size = sub_dim * sub_dim * sizeof(double);
	MPI_Alloc_mem(loc_mat_size, MPI_INFO_NULL, &loc_A);
	MPI_Alloc_mem(loc_mat_size, MPI_INFO_NULL, &loc_B);
	MPI_Alloc_mem(loc_mat_size, MPI_INFO_NULL, &loc_C);
	for (int i = 0; i < sub_dim * sub_dim; ++i) {
		loc_C[i] = 0.0;
	}	

	if (rank == 0) {
		// Read matrices
		/* for (int i = 0; i < m_dim_nopad; ++i) {
			for (int j = 0; j < m_dim_nopad; ++j) {
				scanf("%lf", &(full_A[I_GMAT(i,j)]));
			}
		}
		for (int i = 0; i < m_dim_nopad; ++i) {
			for (int j = 0; j < m_dim_nopad; ++j) {
				scanf("%lf", &(full_B[I_GMAT(i,j)]));
			}
		} */

		for (int i = 0; i < m_dim_nopad; ++i) {
			for (int j = 0; j < m_dim_nopad; ++j) {
				if (i == j) {
					full_A[I_GMAT(i,j)]=1;
					full_B[I_GMAT(i,j)]=1;
				} else {
					full_A[I_GMAT(i,j)]=0;
					full_B[I_GMAT(i,j)]=0;
				}
			}
		}


		start = MPI_Wtime();
	}

	MPI_Datatype submat;
	MPI_Type_vector(sub_dim, sub_dim, m_dim, MPI_DOUBLE, &submat);
	MPI_Type_commit(&submat);

	MPI_Win win_A, win_B;
	MPI_Win_create(full_A, m_dim * m_dim * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_A);
	MPI_Win_create(full_B, m_dim * m_dim * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_B);
	
	for (int shift = 0; shift < p_dim; ++shift) {
		MPI_Win_fence(0, win_A);
		MPI_Get(loc_A, sub_dim * sub_dim, MPI_DOUBLE, 0, my_row * m_dim * sub_dim + ((my_col - my_row - shift + 2*p_dim) % p_dim) * sub_dim, 1, submat, win_A);
		MPI_Win_fence(0, win_A);

		MPI_Win_fence(0, win_B);
		MPI_Get(loc_B, sub_dim * sub_dim, MPI_DOUBLE, 0, ((my_row - my_col - shift + 2*p_dim) % p_dim) * m_dim * sub_dim + my_col * sub_dim, 1, submat, win_B);
		MPI_Win_fence(0, win_B);

		// Matrix multiplication
		for (int i = 0; i < sub_dim; ++i) {
			for (int j = 0; j < sub_dim; ++j) {
				for (int k = 0; k < sub_dim; ++k) {
					loc_C[I_SMAT(i,j)] += loc_A[I_SMAT(i,k)] * loc_B[I_SMAT(k,j)];
				}
			}
		}
	}
	MPI_Win_free(&win_A);
	MPI_Win_free(&win_B);

	MPI_Win win_C;
	MPI_Win_create(full_C, m_dim * m_dim * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_C);
	MPI_Win_fence(0, win_C);
	MPI_Put(loc_C, sub_dim * sub_dim, MPI_DOUBLE, 0, my_row * m_dim * sub_dim + my_col * sub_dim, 1, submat, win_C);
	MPI_Win_fence(0, win_C);
	MPI_Win_free(&win_C);

	MPI_Free_mem(loc_A);
	MPI_Free_mem(loc_B);
	MPI_Free_mem(loc_C);

	if (rank ==0) {
		end = MPI_Wtime();
		printf("Execution time: %f\n", end - start);
	}

	MPI_Free_mem(full_A);
	MPI_Free_mem(full_B);
	MPI_Free_mem(full_C);
	MPI_Finalize();
	return 0;
main_err: {
	     MPI_Finalize();
	     return 1;
     }
}
