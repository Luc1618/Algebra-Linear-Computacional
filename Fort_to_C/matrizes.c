#include <stdio.h> 
#include <stdlib.h>  
#include "matrix_A.h"
#include "vector_b.h"


extern void dgeqrf_ (const int *m, const int *n, double *A, const int *ldA, double *tau, double *work, const int *lwork, int *info);

extern void dormqr_ (const char *side, const char *trans, const int *m, const int *n, const int *k, double *A, const int *ldA, double *tau, double *c, const int *ldC, double *work, int *lwork, int *info );

extern void dtrsv_ (const char *uplo, const char *trans, const char *diag, const int *n, double *a, const int *lda, double *x, const int *incx ) ;

int main(int argc, char **argv)
{
	int m = 20 ;
	int n = 20 ;
	double tau[ldA] ;
	int lwork = -1;
	double work0[2] ;
	double *work ;
	int info ;
	int one = 1 ; 
	int k = n ;
    int i ;
	double x[n];

dgeqrf_(&m, &n, A, &ldA, tau, work0, &lwork, &info) ;

lwork = work0[0] ;

printf("Tamanho ótimo de área de trabalho=%d\n", lwork) ;

work = calloc(lwork, sizeof(double)) ;

/* FATORAÇÃO QR*/
dgeqrf_(&m, &n, A, &ldA, tau, work0, &lwork, &info) ;

printf("Após dgeqrf info=%d\n", info) ;

/* APLICA Q^T a b*/
dormqr_ ("L", "T", &m, &one, &k, A, &ldA, tau, b, &ldb, work, &lwork, &info) ;

printf("Após dormqr info=%d\n", info) ;

/* RESOLVE R x = Qᵀ b */
dtrsv_ ("U", "N", "N", &n, A, &ldA, b, &one) ;

for (i = 0; i < n; i++) /*x[i] = b[i]*/

    printf("%-16e %-16e %-16e\n", x[i], b[i], fabs((x[i] - b[i]))/b[i]) ;

free(work) ;

} ;


