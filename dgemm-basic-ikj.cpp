const char* dgemm_desc = "Basic implementation, ikj-loop dgemm.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < n; j++)
        C[i * n + j] += A[i * n + k] * B[k * n + j];
}
