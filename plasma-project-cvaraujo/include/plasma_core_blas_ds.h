/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/include/plasma_core_blas_zc.h, mixed zc -> ds, Fri Sep 28 17:38:01 2018
 *
 **/
#ifndef PLASMA_CORE_BLAS_DS_H
#define PLASMA_CORE_BLAS_DS_H

#include "plasma_async.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
void plasma_core_dlag2s(int m, int n,
                 double *A,  int lda,
                 float *As, int ldas);

void plasma_core_slag2d(int m, int n,
                 float *As, int ldas,
                 double *A,  int lda);

/******************************************************************************/
void plasma_core_omp_dlag2s(int m, int n,
                     double *A,  int lda,
                     float *As, int ldas,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slag2d(int m, int n,
                     float *As, int ldas,
                     double *A,  int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // PLASMA_CORE_BLAS_DS_H
