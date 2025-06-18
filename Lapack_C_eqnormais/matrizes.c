#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "matrizA.h"
#include "vetorB.h"
#include "vetorX.h"

extern void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *A, const int *ldA, double *beta, double *C, int *ldc);

extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, const int *ldA, double *x, int *incx, double *beta, double *y, int *incy);

extern void dposv_(char *uplo, int *n, int *nrhs, double *A, const int *ldA, double *B, int *ldb, int *info);

int main() {

    // Não esquer de alterar ao mudar de exemplo
    int m = 100;
    int n = 15;
    // Alterar tamanho 
    double AtA[n*n];  // n x n
    double Atb[n*n];  // n x 1
    double alpha = 1.0;
    double beta = 0.0;
    int k = m;
    int l = n ;
    int info;
    int inc1 = 1;
    int i;

printf("m = %d, n = %d, lda = %d\n", m, n, ldA);

// 1. AtA = Aᵗ A
dsyrk_("U", "T", &l, &k, &alpha, A, &ldA, &beta, AtA, &n);

// 2. Atb = Aᵗ b
dgemv_("T", &m, &n, &alpha, A, &ldA, b, &inc1, &beta, Atb, &inc1);

// 3. Resolve AtA x = Atb usando Cholesky
dposv_("U", &n, &inc1, AtA, &n, Atb, &n, &info);

printf("Erro: dposv_ retornou info = %d\n", info);

// Solução Julia | Solução C/lapack | Erro relativo
printf("Solução Julia | Solução C/lapack | Erro relativo \n");
for (i = 0; i < n; i++)

    printf("%-16e %-16e %-16e\n", x[i], Atb[i], fabs((x[i] - Atb[i]))/Atb[i]) ;


}




