#include "dgemm-strassen-winograd.hpp"

const char* dgemm_desc = "Strassen-Winograd dgemm 256.";

void square_dgemm(int n, double* A, double* B, double* C) 
{
  Strassen<256>(n, A, B, C);
}