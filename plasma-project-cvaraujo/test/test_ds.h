/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/test/test_zc.h, mixed zc -> ds, Fri Sep 28 17:38:01 2018
 *
 **/
#ifndef TEST_DS_H
#define TEST_DS_H

#include "test.h"

//==============================================================================
// test routines
//==============================================================================
void test_dsgesv(param_value_t param[], bool run);
void test_dsposv(param_value_t param[], bool run);
void test_dsgbsv(param_value_t param[], bool run);
void test_dlag2s(param_value_t param[], bool run);
void test_slag2d(param_value_t param[], bool run);

#endif // TEST_DS_H
