#include "dgemm-strassen-winograd.hpp"

const char* dgemm_desc = "Strassen-Winograd dgemm 64.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  std::vector<double> buf(n * n);
  _C = buf.data() + 0;
  Strassen<16>(n, A, B, _C);
  matrixAdd(n, C, _C, C);
}