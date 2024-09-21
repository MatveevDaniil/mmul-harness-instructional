#include "dgemm-strassen-winograd.hpp"

const char* dgemm_desc = "Strassen-Winograd dgemm 16.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  Strassen<16>(n, A, B, C);
}