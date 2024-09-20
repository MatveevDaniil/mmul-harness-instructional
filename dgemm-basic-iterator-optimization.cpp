const char* dgemm_desc = "Basic implementation, three-loop dgemm with iterator optimization.";

void square_dgemm(int n, double* A, double* B, double* C) 
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