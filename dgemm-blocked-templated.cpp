#include <vector>
#include <iostream>
#include <stdio.h>

const char* dgemm_desc = "Blocked dgemm.";

template<int n>
void block_dgemm(double* A, double* B, double* C) 
{
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < n; j++)
        C[i * n + j] += A[i * n + k] * B[k * n + j];
}

// void print_matrix(int n, double *A) {
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++)
//       printf("%f ", A[i + j * n]);
//     printf("\n");
//   }
// }

template<int block_size>
void copy_to_block(int n, double *M, double *block, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1)
      *B_ij = *M_ij;
}

template<int block_size>
void add_from_block(int n, double *M, double *block, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1) {
      *M_ij += *B_ij;
      *B_ij = 0;
    }
}

template<int block_size>
void square_dgemm_blocked_templated(int n, double* A, double* B, double* C) 
{
  std::vector<double> buf(3 * block_size * block_size);
  double* A_block = buf.data() + 0;
  double* B_block = A_block + block_size * block_size;
  double* C_block = B_block + block_size * block_size;
  for (int i = 0; i < n / block_size; i++)
    for (int j = 0; j < n / block_size; j++)
      for (int k = 0; k < n / block_size; k++) {
        copy_to_block<block_size>(n, A, A_block, i, k);
        copy_to_block<block_size>(n, B, B_block, k, j);
        block_dgemm<block_size>(A_block, B_block, C_block);
        add_from_block<block_size>(n, C, C_block, i, j);
      }
}

void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
  switch (block_size) {
    case 2:
      square_dgemm_blocked_templated<2>(n, A, B, C);
      break;
    case 16:
      square_dgemm_blocked_templated<16>(n, A, B, C);
      break;
    case 32:
      square_dgemm_blocked_templated<32>(n, A, B, C);
      break;
    case 64:
      square_dgemm_blocked_templated<64>(n, A, B, C);
      break;
    default:
      std::cerr << "Unsupported block size" << std::endl;
  }
}
