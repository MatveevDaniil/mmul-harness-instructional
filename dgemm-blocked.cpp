#include <vector>
#include <iostream>
#include <stdio.h>

const char* dgemm_desc = "Blocked dgemm.";

void block_dgemm(int n, double* A, double* B, double* C) 
{
  double dot_prod, *A_i, *C_i, *B_j, *C_ij, *A_ik, *B_kj;
  for (A_i = A, C_i = C; A_i < A + n * n; A_i += n, C_i += n)
    for (B_j = B, C_ij = C_i; B_j < B + n; B_j += 1, C_ij += 1) {
      dot_prod = 0;
      for (A_ik = A_i, B_kj = B_j; B_kj < B_j + n * n; A_ik += 1, B_kj += n)
        dot_prod += (*A_ik) * (*B_kj);
      *C_ij += dot_prod;
    }
}

void print_matrix(int n, double *A) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%f ", A[i + j * n]);
    printf("\n");
  }
}

void copy_to_block(int n, int block_size, double *M, double *block, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1)
      *B_ij = *M_ij;
}

void add_from_block(int n, int block_size, double *M, double *block, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1) {
      *M_ij += *B_ij;
      *B_ij = 0;
    }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
  std::vector<double> buf(3 * block_size * block_size);
  double* A_block = buf.data() + 0;
  double* B_block = A_block + block_size * block_size;
  double* C_block = B_block + block_size * block_size;
  for (int i = 0; i < n / block_size; i++)
    for (int j = 0; j < n / block_size; j++)
      for (int k = 0; k < n / block_size; k++) {
        copy_to_block(n, block_size, A, A_block, i, k);
        copy_to_block(n, block_size, B, B_block, k, j);
        block_dgemm(block_size, A_block, B_block, C_block);
        add_from_block(n, block_size, C, C_block, i, j);
      }
}
