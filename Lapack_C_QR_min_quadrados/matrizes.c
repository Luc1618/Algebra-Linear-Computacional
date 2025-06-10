#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "matrizA.h"
#include "vetorB.h"
#include "vetorX.h"

extern void dgeqrf_(const int *m, const int *n, double *A, const int *ldA, double *tau, double *work, const int *lwork, int *info);

extern void dormqr_(const char *side, const char *trans, const int *m, const int *n, const int *k, double *A, const int *ldA, double *tau, double *c, const int *ldC, double *work, int *lwork, int *info );

extern void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *n, double *a, const int *lda, double *x, const int *incx ) ;

int main(int argc, char **argv)
{	
	// Não esquecer de alterar quando
	// mudar de exemplo !!!
	int m = 100 ;
	int n = 7 ;
	double tau[ldA] ;
	int lwork = -1;
	double work0[2] ;
	double *work ;
	int info ;
	int one = 1 ; 
	int k = n ;
    	int i ;


// Decomposição Qr - chamada teste (work0)
dgeqrf_ (&m, &n, A, &ldA, tau, work0, &lwork, &info) ;

lwork = work0[0] ;

printf("Tamanho ótimo de área de trabalho=%d\n", lwork) ;

work = calloc(lwork, sizeof(double)) ;

// Decomposição Qr - atualização work 
dgeqrf_ (&m, &n, A, &ldA, tau, work, &lwork, &info) ;

printf("Após dgeqrf info=%d\n", info) ;

// aplica uma Q (mxm) ortogonal da dgeqrf
// Q^T(mxm)b(mx1)
dormqr_ ("L", "T", &m, &one, &k, A, &ldA, tau, b, &ldb, work, &lwork, &info) ;

printf("Após dormqr info=%d\n", info) ;

// soluciona sistema upper triangular
dtrsv_ ("U", "N", "N", &n, A, &ldA, b, &one) ;


// Solução Julia | Solução C/lapack | Erro relativo
printf("Solução Julia | Solução C/lapack | Erro relativo \n");
for (i = 0; i < n; i++)

    printf("%-16e %-16e %-16e\n", x[i], b[i], fabs((x[i] - b[i]))/b[i]) ;

free(work) ;

} ;


