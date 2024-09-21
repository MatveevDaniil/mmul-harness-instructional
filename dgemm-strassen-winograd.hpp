// This code is modification of the following implementation: 
// https://github.com/lrvine/Strassen-algorithm/blob/master/strassen.c

void seq_dgemm(int n, double* A, double* B, double* C) 
{
  for (int k = 0; k < n; k++)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void copy_to_block(int n, double *M, double *block, int block_size, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1)
      *B_ij = *M_ij;
}

void add_from_block(int n, double *M, double *block, int block_size, int block_i, int block_j) {
  M += block_i * block_size * n + block_j * block_size;
  double *M_i, *B_i, *M_ij, *B_ij;
  for (M_i = M, B_i = block; M_i < M + block_size * n; M_i += n, B_i += block_size)
    for (M_ij = M_i, B_ij = B_i; M_ij < M_i + block_size; M_ij += 1, B_ij += 1)
      *M_ij += *B_ij;
}

void matrixAdd(int n, double *A, double *B, double *C) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i * n + j] = A[i * n + j] + B[i * n + j];
}

void matrixSub(int n, double *A, double *B, double *C) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i * n + j] = A[i * n + j] - B[i * n + j];
}

template <int seq_limit>
void Strassen(int n, double* X, double* Y, double* Z) 
{
  if (n <= seq_limit) {
    seq_dgemm(n, X, Y, C);
    return;
  }
  int n_2 = n / 2;
  int n2_sq = n_2 * n_2;

  std::vector<double> a(n2_sq);
  std::vector<double> b(n2_sq);
  std::vector<double> c(n2_sq);
  std::vector<double> d(n2_sq);
  std::vector<double> A(n2_sq);
  std::vector<double> C(n2_sq);
  std::vector<double> B(n2_sq);
  std::vector<double> D(n2_sq);
  std::vector<double> t(n2_sq, 0);
  std::vector<double> u(n2_sq, 0);
  std::vector<double> v(n2_sq, 0);
  std::vector<double> w(n2_sq, 0);

  copy_to_block(n, X, a, n_2, 0, 0);
  copy_to_block(n, X, b, n_2, 0, n_2);
  copy_to_block(n, X, c, n_2, n_2, 0);
  copy_to_block(n, X, d, n_2, n_2, n_2);
  copy_to_block(n, Y, A, n_2, 0, 0);
  copy_to_block(n, Y, C, n_2, 0, n_2);
  copy_to_block(n, Y, B, n_2, n_2, 0);
  copy_to_block(n, Y, D, n_2, n_2, n_2);

  std::vector<double> tempA(n2_sq, 0);
  std::vector<double> tempB(n2_sq, 0);

  matrixSub(n_2, a, c, tempA);
  matrixAdd(n_2, c, d, c);
  matrixSub(n_2, C, A, tempB);
  matrixSub(n_2, D, C, C);
  Strassen<seq_limit>(n_2, tempA, C, v);
  matrixSub(n_2, c, a, tempA);
  Strassen<seq_limit>(n_2, a, A, t);
  matrixSub(n_2, D, tempB, A);
  Strassen<seq_limit>(n_2, c, tempB, a);
  matrixSub(n_2, A, B, tempB);
  Strassen<seq_limit>(n_2, d, tempB, c);

  matrixSub(n_2, b, tempA, d);
  Strassen<seq_limit>(n_2, tempA, A, tempB);
  matrixAdd(n_2, t, tempB, w);
  Strassen<seq_limit>(n_2, b, B, tempA);
  matrixAdd(n_2, t, tempA, t);
  matrixAdd(n_2, w, a, tempA);
  matrixAdd(n_2, w, v, w);
  matrixSub(n_2, w, c, v);
  matrixAdd(n_2, w, a, w);
  Strassen<seq_limit>(n_2, d, D, b);
  matrixAdd(n_2, tempA, b, u);

  four_blocks_to_matrix(n_2, padding, Z, t, u, v, w);
}