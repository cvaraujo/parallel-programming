/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/include/core_lapack_z.h, normal z -> d, Fri Sep 28 17:38:00 2018
 *
 **/

#ifndef PLASMA_CORE_LAPACK_D_H
#define PLASMA_CORE_LAPACK_D_H

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK_GLOBAL is Fortran name mangling macro from LAPACKE

// LAPACKE_dlantr broken (returns 0) in LAPACKE < 3.6.1
#ifndef LAPACK_dlantr
#define LAPACK_dlantr LAPACK_GLOBAL(dlantr, DLANTR)
double LAPACK_dlantr(const char *norm, const char *uplo, const char *diag,
                     const lapack_int *m, const lapack_int *n,
                     const double *A, const lapack_int *lda,
                     double *work);
#endif

// LAPACKE_dlascl not available in LAPACKE < 3.6.0
#ifndef LAPACK_dlascl
#define LAPACK_dlascl LAPACK_GLOBAL(dlascl, DLASCL)
void LAPACK_dlascl(const char *type, const lapack_int *kl, const lapack_int *ku,
                   const double *cfrom, const double *cto,
                   const lapack_int *m, const lapack_int *n,
                   double *A, const lapack_int *lda,
                   lapack_int *info);
#endif

// LAPACKE_dlassq not available yet
#ifndef LAPACK_dlassq
#define LAPACK_dlassq LAPACK_GLOBAL(dlassq, DLASSQ)
void LAPACK_dlassq(const lapack_int *n, const double *x, const lapack_int *incx,
                   double *scale, double *sumsq);
#endif

// LAPACKE_dlangb not available yet
#ifndef LAPACK_dlangb
#define LAPACK_dlangb LAPACK_GLOBAL(dlangb, DLANGB)
double LAPACK_dlangb(const char *norm,
                     const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                     const double *A, const lapack_int *lda,
                     double *work);

#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // PLASMA_CORE_LAPACK_D_H
