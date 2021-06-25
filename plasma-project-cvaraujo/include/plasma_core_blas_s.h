/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/include/plasma_core_blas_z.h, normal z -> s, Fri Sep 28 17:38:01 2018
 *
 **/
#ifndef PLASMA_CORE_BLAS_S_H
#define PLASMA_CORE_BLAS_S_H

#include "plasma_async.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "plasma_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REAL

/******************************************************************************/
#ifdef COMPLEX
float fabsf(float alpha);
#endif

int plasma_core_sgeadd(plasma_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

int plasma_core_sgelqt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void plasma_core_sgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

int plasma_core_sgeqrt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void plasma_core_sgessq(int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void plasma_core_sgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                 volatile int *max_idx, volatile float *max_val,
                 volatile int *info, plasma_barrier_t *barrier);

int plasma_core_ssygst(int itype, plasma_enum_t uplo,
                int n,
                float *A, int lda,
                float *B, int ldb);

void plasma_core_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void plasma_core_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 float alpha, const float *A, int lda,
                                           const float *B, int ldb,
                 float beta,                    float *C, int ldc);

void plasma_core_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

void plasma_core_ssyssq(plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void plasma_core_ssyssq(plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void plasma_core_slacpy(plasma_enum_t uplo, plasma_enum_t transa,
                 int m, int n,
                 const float *A, int lda,
                       float *B, int ldb);

void plasma_core_slacpy_lapack2tile_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *A, int lda,
                                        float *B, int ldb);

void plasma_core_slacpy_tile2lapack_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *B, int ldb,
                                        float *A, int lda);

void plasma_core_slange(plasma_enum_t norm,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *result);

void plasma_core_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void plasma_core_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void plasma_core_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *value);

void plasma_core_slascl(plasma_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 float *A, int lda);

void plasma_core_slaset(plasma_enum_t uplo,
                 int m, int n,
                 float alpha, float beta,
                 float *A, int lda);

void plasma_core_sgeswp(plasma_enum_t colrow,
                 plasma_desc_t A, int k1, int k2, const int *ipiv, int incx);

void plasma_core_ssyswp(int rank, int num_threads,
                 int uplo, plasma_desc_t A, int k1, int k2, const int *ipiv,
                 int incx, plasma_barrier_t *barrier);

int plasma_core_slauum(plasma_enum_t uplo,
                int n,
                float *A, int lda);

int plasma_core_spamm(plasma_enum_t op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const float *A1, int lda1,
                     float *A2, int lda2,
               const float *V,  int ldv,
                     float *W,  int ldw);

int plasma_core_sparfb(plasma_enum_t side, plasma_enum_t trans, plasma_enum_t direct,
                plasma_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int plasma_core_spemv(plasma_enum_t trans, int storev,
               int m, int n, int l,
               float alpha,
               const float *A, int lda,
               const float *X, int incx,
               float beta,
               float *Y, int incy,
               float *work);

int plasma_core_spotrf(plasma_enum_t uplo,
                int n,
                float *A, int lda);

void plasma_core_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void plasma_core_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void plasma_core_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

int plasma_core_stradd(plasma_enum_t uplo, plasma_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

void plasma_core_strmm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void plasma_core_strsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void plasma_core_strssq(plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

int plasma_core_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                int n,
                float *A, int lda);

int plasma_core_stslqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int plasma_core_stsmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int plasma_core_stsmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int plasma_core_stsqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int plasma_core_sttlqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int plasma_core_sttmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int plasma_core_sttmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int plasma_core_sttqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int plasma_core_sormlq(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

int plasma_core_sormqr(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

/******************************************************************************/
void plasma_core_omp_samax(int colrow, int m, int n,
                     const float *A, int lda,
                     float *values,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgeadd(
    plasma_enum_t transa, int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgelqt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgeqrt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgessq(int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_ssygst(int itype, plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     float *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,                    float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                    int n, int k,
                    float alpha, const float *A, int lda,
                    float beta,        float *C, int ldc,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyssq(plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyssq(plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_slacpy(plasma_enum_t uplo, plasma_enum_t transa,
                     int m, int n,
                     const float *A, int lda,
                           float *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slacpy_lapack2tile_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *A, int lda,
                                            float *B, int ldb);

void plasma_core_omp_slacpy_tile2lapack_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *B, int ldb,
                                            float *A, int lda);

void plasma_core_omp_slange(plasma_enum_t norm,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *result,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slange_aux(plasma_enum_t norm,
                         int m, int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slantr_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         plasma_enum_t diag,
                         int m, int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_slascl(plasma_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     float *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_slaset(plasma_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     float alpha, float beta,
                     float *A);

void plasma_core_omp_slauum(plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_spotrf(plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_ssyrk(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_stradd(
    plasma_enum_t uplo, plasma_enum_t transa,
    int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_strmm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_strsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_strssq(plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                     int n,
                     float *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_stslqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_stsmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_stsmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_stsqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sttlqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sttmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sttmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sttqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sormlq(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_sormqr(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#undef REAL

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // PLASMA_CORE_BLAS_S_H
