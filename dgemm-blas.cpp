#include <cblas.h>

const char* dgemm_desc = "Reference dgemm.";
void square_dgemm(int n, double* A, double* B, double* C) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., A, n, B, n, 1., C, n);
}
