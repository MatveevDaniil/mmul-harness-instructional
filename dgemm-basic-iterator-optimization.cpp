const char* dgemm_desc = "Basic implementation, three-loop dgemm with iterator optimization.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  double *A_i, *C_i, *B_k, *C_ij, *A_ik, *B_kj;
  for (A_i = A, C_i = C; A_i < A + n * n; A_i += n, C_i += n)
    for (B_k = B, A_ik = A_i; A_ik < A_i + n;  A_ik++, B_k += n)
      for (C_ij = C_i, B_kj = B_k; C_ij < C_i + n; C_ij++, B_kj++)
        *C_ij += (*A_ik) * (*B_kj);
}