/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/include/plasma_core_blas_z.h, normal z -> d, Fri Sep 28 17:38:01 2018
 *
 **/
#ifndef PLASMA_CORE_BLAS_D_H
#define PLASMA_CORE_BLAS_D_H

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
double fabs(double alpha);
#endif

int plasma_core_dgeadd(plasma_enum_t transa,
                int m, int n,
                double alpha, const double *A, int lda,
                double beta,        double *B, int ldb);

int plasma_core_dgelqt(int m, int n, int ib,
                double *A, int lda,
                double *T, int ldt,
                double *tau,
                double *work);

void plasma_core_dgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

int plasma_core_dgeqrt(int m, int n, int ib,
                double *A, int lda,
                double *T, int ldt,
                double *tau,
                double *work);

void plasma_core_dgessq(int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

void plasma_core_dgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                 volatile int *max_idx, volatile double *max_val,
                 volatile int *info, plasma_barrier_t *barrier);

int plasma_core_dsygst(int itype, plasma_enum_t uplo,
                int n,
                double *A, int lda,
                double *B, int ldb);

void plasma_core_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

void plasma_core_dsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 double alpha, const double *A, int lda,
                                           const double *B, int ldb,
                 double beta,                    double *C, int ldc);

void plasma_core_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                double alpha, const double *A, int lda,
                double beta,        double *C, int ldc);

void plasma_core_dsyssq(plasma_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

void plasma_core_dsyssq(plasma_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

void plasma_core_dlacpy(plasma_enum_t uplo, plasma_enum_t transa,
                 int m, int n,
                 const double *A, int lda,
                       double *B, int ldb);

void plasma_core_dlacpy_lapack2tile_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const double *A, int lda,
                                        double *B, int ldb);

void plasma_core_dlacpy_tile2lapack_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const double *B, int ldb,
                                        double *A, int lda);

void plasma_core_dlange(plasma_enum_t norm,
                 int m, int n,
                 const double *A, int lda,
                 double *work, double *result);

void plasma_core_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *work, double *value);

void plasma_core_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *work, double *value);

void plasma_core_dlantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const double *A, int lda,
                 double *work, double *value);

void plasma_core_dlascl(plasma_enum_t uplo,
                 double cfrom, double cto,
                 int m, int n,
                 double *A, int lda);

void plasma_core_dlaset(plasma_enum_t uplo,
                 int m, int n,
                 double alpha, double beta,
                 double *A, int lda);

void plasma_core_dgeswp(plasma_enum_t colrow,
                 plasma_desc_t A, int k1, int k2, const int *ipiv, int incx);

void plasma_core_dsyswp(int rank, int num_threads,
                 int uplo, plasma_desc_t A, int k1, int k2, const int *ipiv,
                 int incx, plasma_barrier_t *barrier);

int plasma_core_dlauum(plasma_enum_t uplo,
                int n,
                double *A, int lda);

int plasma_core_dpamm(plasma_enum_t op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const double *A1, int lda1,
                     double *A2, int lda2,
               const double *V,  int ldv,
                     double *W,  int ldw);

int plasma_core_dparfb(plasma_enum_t side, plasma_enum_t trans, plasma_enum_t direct,
                plasma_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int plasma_core_dpemv(plasma_enum_t trans, int storev,
               int m, int n, int l,
               double alpha,
               const double *A, int lda,
               const double *X, int incx,
               double beta,
               double *Y, int incy,
               double *work);

int plasma_core_dpotrf(plasma_enum_t uplo,
                int n,
                double *A, int lda);

void plasma_core_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

void plasma_core_dsyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void plasma_core_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                double alpha, const double *A, int lda,
                double beta,        double *C, int ldc);

int plasma_core_dtradd(plasma_enum_t uplo, plasma_enum_t transa,
                int m, int n,
                double alpha, const double *A, int lda,
                double beta,        double *B, int ldb);

void plasma_core_dtrmm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                double alpha, const double *A, int lda,
                                                double *B, int ldb);

void plasma_core_dtrsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                double alpha, const double *A, int lda,
                                                double *B, int ldb);

void plasma_core_dtrssq(plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

int plasma_core_dtrtri(plasma_enum_t uplo, plasma_enum_t diag,
                int n,
                double *A, int lda);

int plasma_core_dtslqt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int plasma_core_dtsmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int plasma_core_dtsmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int plasma_core_dtsqrt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int plasma_core_dttlqt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int plasma_core_dttmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int plasma_core_dttmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int plasma_core_dttqrt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int plasma_core_dormlq(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const double *A,    int lda,
                const double *T,    int ldt,
                      double *C,    int ldc,
                      double *work, int ldwork);

int plasma_core_dormqr(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const double *A,    int lda,
                const double *T,    int ldt,
                      double *C,    int ldc,
                      double *work, int ldwork);

/******************************************************************************/
void plasma_core_omp_damax(int colrow, int m, int n,
                     const double *A, int lda,
                     double *values,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgeadd(
    plasma_enum_t transa, int m, int n,
    double alpha, const double *A, int lda,
    double beta,        double *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgelqt(int m, int n, int ib,
                     double *A, int lda,
                     double *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgeqrt(int m, int n, int ib,
                     double *A, int lda,
                     double *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgessq(int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dgessq_aux(int n,
                         const double *scale, const double *sumsq,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dsygst(int itype, plasma_enum_t uplo,
                     int n,
                     double *A, int lda,
                     double *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,                    double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                    int n, int k,
                    double alpha, const double *A, int lda,
                    double beta,        double *C, int ldc,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyssq(plasma_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *scale, double *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyssq(plasma_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *scale, double *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyssq_aux(int m, int n,
                         const double *scale, const double *sumsq,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dlacpy(plasma_enum_t uplo, plasma_enum_t transa,
                     int m, int n,
                     const double *A, int lda,
                           double *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlacpy_lapack2tile_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const double *A, int lda,
                                            double *B, int ldb);

void plasma_core_omp_dlacpy_tile2lapack_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const double *B, int ldb,
                                            double *A, int lda);

void plasma_core_omp_dlange(plasma_enum_t norm,
                     int m, int n,
                     const double *A, int lda,
                     double *work, double *result,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlange_aux(plasma_enum_t norm,
                         int m, int n,
                         const double *A, int lda,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *work, double *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const double *A, int lda,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *work, double *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const double *A, int lda,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dlantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const double *A, int lda,
                     double *work, double *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlantr_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         plasma_enum_t diag,
                         int m, int n,
                         const double *A, int lda,
                         double *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_core_omp_dlascl(plasma_enum_t uplo,
                     double cfrom, double cto,
                     int m, int n,
                     double *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dlaset(plasma_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     double alpha, double beta,
                     double *A);

void plasma_core_omp_dlauum(plasma_enum_t uplo,
                     int n,
                     double *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dpotrf(plasma_enum_t uplo,
                     int n,
                     double *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dsyrk(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
    double beta,        double *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtradd(
    plasma_enum_t uplo, plasma_enum_t transa,
    int m, int n,
    double alpha, const double *A, int lda,
    double beta,        double *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtrmm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    double alpha, const double *A, int lda,
                                    double *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtrsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    double alpha, const double *A, int lda,
                                    double *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtrssq(plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtrtri(plasma_enum_t uplo, plasma_enum_t diag,
                     int n,
                     double *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtslqt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtsmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V,  int ldv,
                     const double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtsmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V, int ldv,
                     const double *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dtsqrt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dttlqt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dttmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V,  int ldv,
                     const double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dttmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V, int ldv,
                     const double *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dttqrt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dormlq(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const double *A, int lda,
                     const double *T, int ldt,
                           double *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_core_omp_dormqr(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const double *A, int lda,
                     const double *T, int ldt,
                           double *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#undef REAL

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // PLASMA_CORE_BLAS_D_H
