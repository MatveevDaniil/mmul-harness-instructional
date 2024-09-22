const char* dgemm_desc = "Basic implementation, ijk-loop dgemm.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        C[i * n + j] += A[i * n + k] * B[k * n + j];
}
