/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/core_blas/core_zlascl.c, normal z -> c, Fri Sep 28 17:38:21 2018
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void plasma_core_clascl(plasma_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 plasma_complex32_t *A, int lda)
{
    // LAPACKE_clascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_clascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void plasma_core_omp_clascl(plasma_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     plasma_complex32_t *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(inout:A[0:lda*n])
    {
        if (sequence->status == PlasmaSuccess)
            plasma_core_clascl(uplo,
                        cfrom, cto,
                        m, n,
                        A, lda);
    }
}
