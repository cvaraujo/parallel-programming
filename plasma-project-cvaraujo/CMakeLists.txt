cmake_minimum_required( VERSION 3.0 FATAL_ERROR )

project( PLASMA VERSION 19.8.1 LANGUAGES C )

set(CMAKE_SUPPRESS_REGENERATION on)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#set( CMAKE_THREAD_PREFER_PTHREAD 1 )
#find_package( Threads )

macro( use_c99 )
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ( CMAKE_C_COMPILER_ID )
  else ( )
    set ( CMAKE_C_STANDARD 99 )
  endif ( )
endmacro( use_c99 )

# PLASMA uses C99 features (in-loop definition of for loop variables)
use_c99()

# use standard module to find LAPACK
find_package( BLAS REQUIRED ) # use BLA_VENDOR (by default "All") "Intel" "Intel10_32" "Intel10_64lp" "Intel10_64p_seq" "IBMESSL" "Generic"

# use standard module to find LAPACK (this seems to call FindBLAS)
find_package( LAPACK REQUIRED )

message( STATUS "BLAS(${BLAS_FOUND})" )  # linker flags `${BLA_LINKER_FLAGS}'" )
message( STATUS "LAPACK(${LAPACK_FOUND})" )  # libraries: ${LAPACK_LIBRARIES}" )

if (BLAS_FOUND)
  # BLA_VENDOR is the input variable to allow only checks for specific vendor BLAS
  #if ( BLA_VENDOR STREQUAL "Intel" OR BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "Intel10_64p_seq" )
  #
  if ( BLAS_LIBRARIES MATCHES "mkl_core" )

    message( STATUS "Found Intel MKL" )
    add_definitions( -DPLASMA_WITH_MKL ) # this is command line only
    set( HAVE_MKL 1 )

  elseif ( BLAS_LIBRARIES MATCHES "essl" )
    message( STATUS "Found IBM ESSL" )

    set( HAVE_ESSL 1 )
    set( PLASMA_WITH_ESSL 1 )

  elseif ( BLAS_LIBRARIES MATCHES "openblas" )
    message( STATUS "Found OpenBLAS" )

    add_definitions( -DPLASMA_WITH_OPENBLAS )
    set( PLASMA_WITH_OPENBLAS 1 )

  elseif ( BLAS_LIBRARIES MATCHES "Accelerate" )
    message( STATUS "Found Apple Accelerate Framework" )

    add_definitions( -DPLASMA_WITH_ACCELERATE )
    set( PLASMA_WITH_ACCELERATE 1 )

  else ()
    message( STATUS "Found Generic BLAS" )
    message( STATUS "vendor ${BLA_VENDOR}, linker flags `${BLA_LINKER_FLAGS}'" )
  endif ()

endif (BLAS_FOUND)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/Modules")

cmake_policy(PUSH)
cmake_policy(SET CMP0074 NEW) # allows to use CBLAS_ROOT and LAPACKE_ROOT

find_package( CBLAS REQUIRED )
include_directories(${CBLAS_INCLUDE_DIRS})
if (CBLAS_PROVIDER STREQUAL "mkl")
  add_definitions(-DPLASMA_WITH_MKL)
  set(PLASMA_WITH_MKL TRUE)
elseif (CBLAS_PROVIDER STREQUAL "netlib")
  add_definitions(-DPLASMA_WITH_NETLIB)
  set(PLASMA_WITH_NETLIB TRUE)
endif()

find_package( LAPACKE REQUIRED )
include_directories(${LAPACKE_INCLUDE_DIRS})

cmake_policy(POP)

set(PLASMA_LINALG_LIBRARIES ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES} ${BLAS_LIBRARIES})

if (PLASMA_DETECT_LUA)
  find_package( Lua )

  if ( LUA_FOUND )
    include_directories( ${LUA_INCLUDE_DIR} )
    add_definitions( -DPLASMA_USE_LUA ) # this is command line only
    set( PLASMA_USE_LUA 1 ) # this will be substituted in the config file
  endif()
endif()

# use standard module to find OpenMP
find_package( OpenMP )
if (OPENMP_FOUND)
  #if (OpenMP_C_VERSION  VERSION_LESS 4.5)
  #  message(FATAL_ERROR "OpenMP C version ${OpenMP_C_VERSION} is too old. Version 4.5 or later is required.")
  #endif()
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  # set C++ flags in case C++ compiler is used to compiler PLASMA
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP not found.")
endif()

# Ensure we are using clang
if(NOT "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    message(WARNING "OmpCluster is not compatible with `${CMAKE_C_COMPILER_ID}`!")
endif()

# Remove annoying warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-return-type -Wno-absolute-value")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type -Wno-absolute-value")

option(USE_TAGGING "Build with OmpTracing tagging function" OFF)
if(USE_TAGGING)
  find_package(omptracing CONFIG REQUIRED)
  include_directories(${omptracing_INCLUDE_DIRS})
endif()

option(USE_OMPC "Build for OmpCluster" ON)
if(USE_OMPC)
  if(OPENMP_FOUND)
    if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
      message(STATUS "Building with OmpCluster support")
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp-targets=x86_64-pc-linux-gnu")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=x86_64-pc-linux-gnu")
      add_definitions("-DUSE_OFFLOAD")
    else()
      message("Not building with OmpCluster support: Require to compile with Clang")
    endif()
  else()
    message("Not building with OmpCluster support: OpenMP not found")
  endif()
else()
  message("Not building with OmpCluster support: Manually disabled")
endif()

set(COMPUTE_SOURCE
  compute/clag2z.c compute/dzamax.c compute/scamax.c compute/samax.c
  compute/damax.c compute/pclag2z.c compute/pdzamax.c
  compute/pzdesc2ge.c compute/pzdesc2pb.c compute/pzdesc2tr.c compute/pzgbtrf.c
  compute/pzge2desc.c compute/pzgeadd.c compute/pzgelqf.c compute/pzgelqf_tree.c
  compute/pzgemm.c compute/pzgeqrf.c compute/pzgeqrf_tree.c compute/pzgeswp.c
  compute/pzgetri_aux.c compute/pzhemm.c compute/pzher2k.c
  compute/pzherk.c compute/pzhetrf_aasen.c compute/pzlacpy.c compute/pzlag2c.c
  compute/pzlangb.c compute/pzlange.c compute/pzlanhe.c compute/pzlansy.c
  compute/pzlantr.c compute/pzlascl.c compute/pzlaset.c compute/pzlauum.c
  compute/pzpb2desc.c compute/pzpbtrf.c compute/pzpotrf.c compute/pzsymm.c
  compute/pzsyr2k.c compute/pzsyrk.c compute/pztbsm.c compute/pztr2desc.c
  compute/pztradd.c compute/pztrmm.c compute/pztrsm.c compute/pztrtri.c
  compute/pzunglq.c compute/pzunglq_tree.c compute/pzungqr.c
  compute/pzungqr_tree.c compute/pzunmlq.c compute/pzunmlq_tree.c
  compute/pzunmqr.c compute/pzunmqr_tree.c compute/zcgbsv.c
  compute/zcposv.c compute/zdesc2ge.c compute/zdesc2pb.c compute/zdesc2tr.c
  compute/zgbsv.c compute/zgbtrf.c compute/zgbtrs.c compute/zge2desc.c
  compute/zgeadd.c compute/zgelqf.c compute/zgelqs.c
  compute/zgels.c compute/zgemm.c compute/zgeqrf.c compute/zgeqrs.c
  compute/zgeswp.c compute/zgetri_aux.c
  compute/zgetri.c compute/zgetrs.c compute/zhemm.c compute/zher2k.c
  compute/zherk.c compute/zhesv.c compute/zhetrf.c compute/zhetrs.c
  compute/zlacpy.c compute/clag2z.c compute/zlag2c.c compute/zlangb.c
  compute/zlange.c
  compute/zlanhe.c compute/zlansy.c compute/zlantr.c compute/zlascl.c
  compute/zlaset.c compute/zlauum.c compute/zpb2desc.c compute/zpbsv.c
  compute/zpbtrf.c compute/zpbtrs.c compute/zpoinv.c compute/zposv.c
  compute/zpotrf.c compute/zpotri.c compute/zpotrs.c compute/zsymm.c
  compute/zsyr2k.c compute/zsyrk.c compute/ztr2desc.c compute/ztradd.c
  compute/ztrmm.c compute/ztrsm.c compute/ztrtri.c compute/zunglq.c
  compute/zungqr.c compute/zunmlq.c compute/zunmqr.c compute/cgelqf.c
  compute/cgemm.c compute/cgeqrf.c compute/cpotrf.c compute/cpotrs.c
  compute/csymm.c compute/csyr2k.c compute/csyrk.c compute/ctradd.c
  compute/ctrmm.c compute/ctrsm.c compute/ctrtri.c compute/cunglq.c
  compute/cungqr.c compute/cunmlq.c compute/cunmqr.c compute/dgelqf.c
  compute/dgemm.c compute/dgeqrf.c compute/dorglq.c compute/dorgqr.c
  compute/dormlq.c compute/dormqr.c compute/dpotrf.c compute/dpotrs.c
  compute/dsymm.c compute/dsyr2k.c compute/dsyrk.c compute/dtradd.c
  compute/dtrmm.c compute/dtrsm.c compute/dtrtri.c compute/sgelqf.c
  compute/sgemm.c compute/sgeqrf.c compute/sorglq.c compute/sorgqr.c
  compute/sormlq.c compute/sormqr.c compute/spotrf.c compute/spotrs.c
  compute/ssymm.c compute/ssyr2k.c compute/ssyrk.c compute/stradd.c
  compute/strmm.c compute/strsm.c compute/strtri.c
  compute/dsposv.c compute/dgbsv.c compute/cgbsv.c compute/sgbsv.c
  compute/dgbtrf.c compute/dgbtrs.c compute/cgbtrf.c compute/cgbtrs.c
  compute/sgbtrf.c compute/sgbtrs.c compute/dgeadd.c compute/cgeadd.c
  compute/sgeadd.c
  compute/dgelqs.c compute/cgelqs.c compute/sgelqs.c compute/dgels.c
  compute/cgels.c compute/sgels.c compute/dgeqrs.c compute/cgeqrs.c
  compute/sgeqrs.c compute/dsgbsv.c
  compute/dgetri.c compute/cgetri.c compute/sgetri.c
  compute/dgetri_aux.c compute/cgetri_aux.c compute/sgetri_aux.c
  compute/dgetrs.c compute/cgetrs.c
  compute/sgetrs.c compute/chemm.c compute/cher2k.c
  compute/cherk.c compute/dsytrf.c compute/dsytrs.c compute/chetrf.c
  compute/chetrs.c compute/ssytrf.c compute/ssytrs.c compute/dsysv.c
  compute/chesv.c compute/ssysv.c compute/dlacpy.c compute/clacpy.c
  compute/slacpy.c compute/dlag2s.c compute/slag2d.c compute/dlange.c
  compute/clange.c compute/slange.c compute/clanhe.c compute/dlansy.c
  compute/clansy.c compute/slansy.c compute/dlantr.c compute/clantr.c
  compute/slantr.c compute/dlascl.c compute/clascl.c compute/slascl.c
  compute/dlaset.c compute/claset.c compute/slaset.c compute/dgeswp.c
  compute/cgeswp.c compute/sgeswp.c compute/dlauum.c compute/clauum.c
  compute/slauum.c compute/dpbsv.c compute/cpbsv.c compute/spbsv.c
  compute/dpbtrf.c compute/dpbtrs.c compute/cpbtrf.c compute/cpbtrs.c
  compute/spbtrf.c compute/spbtrs.c compute/dlangb.c compute/clangb.c
  compute/slangb.c compute/dposv.c compute/cposv.c compute/sposv.c
  compute/dpoinv.c compute/cpoinv.c compute/spoinv.c compute/dpotri.c
  compute/cpotri.c compute/spotri.c
  compute/pslange.c compute/pclaset.c compute/psorglq_tree.c
  compute/psormqr_tree.c compute/pdgelqf_tree.c compute/pslag2d.c
  compute/pcunmqr_tree.c compute/psgeqrf_tree.c compute/pspotrf.c
  compute/pdsytrf_aasen.c compute/pslauum.c compute/pssytrf_aasen.c
  compute/pstrsm.c compute/psgeqrf.c compute/pcgelqf_tree.c
  compute/pcunglq_tree.c compute/pctrmm.c compute/pstrtri.c
  compute/pcungqr_tree.c compute/pcsymm.c compute/psormqr.c compute/pdgemm.c
  compute/pdlacpy.c compute/psgeadd.c compute/pdtrmm.c compute/pcungqr.c
  compute/pcgemm.c compute/pslansy.c compute/pdtradd.c compute/pdormqr_tree.c
  compute/pdtbsm.c compute/psormlq.c compute/pdpotrf.c compute/pcunglq.c
  compute/pchemm.c compute/psgeswp.c compute/pcher2k.c compute/pdgetri_aux.c
  compute/pcgeqrf_tree.c compute/pdorglq.c compute/pdlange.c
  compute/pcunmlq_tree.c compute/pdgeqrf.c compute/pdlauum.c
  compute/pdlaset.c compute/pclascl.c compute/pclauum.c compute/pcgeadd.c
  compute/pdorglq_tree.c compute/pdtrsm.c compute/psorglq.c
  compute/pslangb.c compute/pdormlq_tree.c compute/pcherk.c compute/pcpbtrf.c
  compute/psgemm.c compute/pdgeqrf_tree.c compute/pdlascl.c compute/pdsyr2k.c
  compute/pdlantr.c compute/pdgeadd.c compute/pclansy.c compute/psgetri_aux.c
  compute/pclantr.c compute/pstradd.c compute/pcgbtrf.c compute/pcsyrk.c
  compute/pctradd.c compute/psgelqf_tree.c compute/pslantr.c compute/pdlag2s.c
  compute/pslag2d.c
  compute/pchetrf_aasen.c compute/pssymm.c compute/pcunmqr.c compute/pclacpy.c
  compute/pdsyrk.c compute/pcsyr2k.c compute/pdgelqf.c compute/pdamax.c
  compute/pslacpy.c compute/pdormqr.c compute/pctrsm.c compute/pclangb.c
  compute/pdlangb.c compute/pscamax.c compute/pdpbtrf.c compute/pcgeqrf.c
  compute/pdgbtrf.c compute/psamax.c compute/pslascl.c compute/psgbtrf.c
  compute/pdgeswp.c compute/pspbtrf.c compute/pctbsm.c compute/pdorgqr.c
  compute/pcgelqf.c compute/pcpotrf.c compute/pstbsm.c compute/pstrmm.c
  compute/pssyr2k.c compute/pclange.c compute/psorgqr.c compute/psormlq_tree.c
  compute/pssyrk.c compute/pdorgqr_tree.c compute/pdsymm.c compute/pslaset.c
  compute/pdlansy.c compute/pcgeswp.c compute/psorgqr_tree.c compute/pctrtri.c
  compute/pcgetri_aux.c compute/pdormlq.c compute/pcunmlq.c
  compute/pclanhe.c compute/pdtrtri.c compute/psgelqf.c
  compute/zdesc2ge.c compute/zdesc2pb.c compute/zdesc2tr.c
  compute/cdesc2ge.c compute/cdesc2pb.c compute/cdesc2tr.c
  compute/ddesc2ge.c compute/ddesc2pb.c compute/ddesc2tr.c
  compute/sdesc2ge.c compute/sdesc2pb.c compute/sdesc2tr.c
  compute/pzdesc2ge.c compute/pzdesc2pb.c compute/pzdesc2tr.c
  compute/pcdesc2ge.c compute/pcdesc2pb.c compute/pcdesc2tr.c
  compute/pddesc2ge.c compute/pddesc2pb.c compute/pddesc2tr.c
  compute/psdesc2ge.c compute/psdesc2pb.c compute/psdesc2tr.c
  compute/zge2desc.c compute/zpb2desc.c compute/ztr2desc.c
  compute/cge2desc.c compute/cpb2desc.c compute/ctr2desc.c
  compute/dge2desc.c compute/dpb2desc.c compute/dtr2desc.c
  compute/sge2desc.c compute/spb2desc.c compute/str2desc.c
  compute/pzge2desc.c compute/pzpb2desc.c compute/pztr2desc.c
  compute/pcge2desc.c compute/pcpb2desc.c compute/pctr2desc.c
  compute/pdge2desc.c compute/pdpb2desc.c compute/pdtr2desc.c
  compute/psge2desc.c compute/pspb2desc.c compute/pstr2desc.c
  control/constants.c control/context.c control/descriptor.c
)

set(CORE_SOURCE
  core_blas/core_clag2z.c core_blas/core_dcabs1.c core_blas/core_scabs1.c
  core_blas/core_dzamax.c core_blas/core_zgeadd.c core_blas/core_zgelqt.c
  core_blas/core_zgemm.c core_blas/core_zgeqrt.c core_blas/core_zgessq.c
  core_blas/core_zgeswp.c core_blas/core_zgetrf.c
  core_blas/core_zhegst.c core_blas/core_zhemm.c core_blas/core_zher2k.c
  core_blas/core_zherk.c core_blas/core_zhessq.c
  core_blas/core_zheswp.c core_blas/core_zlacpy_band.c core_blas/core_zlacpy.c
  core_blas/core_zlag2c.c core_blas/core_zlange.c
  core_blas/core_zlanhe.c core_blas/core_zlansy.c core_blas/core_zlantr.c
  core_blas/core_zlascl.c core_blas/core_zlaset.c
  core_blas/core_zlauum.c core_blas/core_zpamm.c core_blas/core_zpemv.c
  core_blas/core_zparfb.c core_blas/core_zpemv.c core_blas/core_zpotrf.c
  core_blas/core_zsymm.c core_blas/core_zsyr2k.c core_blas/core_zsyrk.c
  core_blas/core_zsyssq.c core_blas/core_ztradd.c
  core_blas/core_ztrmm.c core_blas/core_ztrsm.c core_blas/core_ztrssq.c
  core_blas/core_ztrtri.c core_blas/core_ztslqt.c
  core_blas/core_ztsmlq.c core_blas/core_ztsmqr.c core_blas/core_ztsqrt.c
  core_blas/core_zttlqt.c core_blas/core_zttmlq.c
  core_blas/core_zttmqr.c core_blas/core_zttqrt.c core_blas/core_zunmlq.c
  core_blas/core_zunmqr.c
  core_blas/core_cgeadd.c core_blas/core_cgemm.c core_blas/core_cgeswp.c
  core_blas/core_cgetrf.c core_blas/core_cheswp.c core_blas/core_clacpy.c
  core_blas/core_clacpy_band.c core_blas/core_cparfb.c core_blas/core_ctrsm.c
  core_blas/core_dgeadd.c core_blas/core_dgemm.c core_blas/core_dgeswp.c
  core_blas/core_dgetrf.c core_blas/core_dlacpy.c core_blas/core_dlacpy_band.c
  core_blas/core_dparfb.c core_blas/core_dsyswp.c core_blas/core_dtrsm.c
  core_blas/core_sgeadd.c core_blas/core_sgemm.c core_blas/core_sgeswp.c
  core_blas/core_sgetrf.c core_blas/core_slacpy.c core_blas/core_slacpy_band.c
  core_blas/core_sparfb.c core_blas/core_ssyswp.c core_blas/core_strsm.c
  core_blas/core_cgelqt.c core_blas/core_cgeqrt.c core_blas/core_cgessq.c
  core_blas/core_chegst.c core_blas/core_chemm.c core_blas/core_cher2k.c
  core_blas/core_cherk.c core_blas/core_chessq.c core_blas/core_clange.c
  core_blas/core_clanhe.c core_blas/core_clansy.c core_blas/core_clantr.c
  core_blas/core_clascl.c core_blas/core_claset.c core_blas/core_clauum.c
  core_blas/core_cpamm.c core_blas/core_cpemv.c core_blas/core_cpotrf.c
  core_blas/core_csymm.c core_blas/core_csyr2k.c core_blas/core_csyrk.c
  core_blas/core_csyssq.c core_blas/core_ctradd.c core_blas/core_ctrmm.c
  core_blas/core_ctrssq.c core_blas/core_ctrtri.c core_blas/core_ctslqt.c
  core_blas/core_ctsmlq.c core_blas/core_ctsmqr.c core_blas/core_ctsqrt.c
  core_blas/core_cttlqt.c core_blas/core_cttmlq.c core_blas/core_cttmqr.c
  core_blas/core_cttqrt.c core_blas/core_cunmlq.c core_blas/core_cunmqr.c
  core_blas/core_damax.c core_blas/core_dgelqt.c core_blas/core_dgeqrt.c
  core_blas/core_dgessq.c core_blas/core_dlag2s.c core_blas/core_dlange.c
  core_blas/core_dlansy.c core_blas/core_dlantr.c core_blas/core_dlascl.c
  core_blas/core_dlaset.c core_blas/core_dlauum.c core_blas/core_dormlq.c
  core_blas/core_dormqr.c core_blas/core_dpamm.c core_blas/core_dpemv.c
  core_blas/core_dpotrf.c core_blas/core_dsygst.c core_blas/core_dsymm.c
  core_blas/core_dsyr2k.c core_blas/core_dsyrk.c core_blas/core_dsyssq.c
  core_blas/core_dtradd.c core_blas/core_dtrmm.c core_blas/core_dtrssq.c
  core_blas/core_dtrtri.c core_blas/core_dtslqt.c core_blas/core_dtsmlq.c
  core_blas/core_dtsmqr.c core_blas/core_dtsqrt.c core_blas/core_dttlqt.c
  core_blas/core_dttmlq.c core_blas/core_dttmqr.c core_blas/core_dttqrt.c
  core_blas/core_samax.c core_blas/core_scamax.c core_blas/core_sgelqt.c
  core_blas/core_sgeqrt.c core_blas/core_sgessq.c core_blas/core_slag2d.c
  core_blas/core_slange.c core_blas/core_slansy.c core_blas/core_slantr.c
  core_blas/core_slascl.c core_blas/core_slaset.c core_blas/core_slauum.c
  core_blas/core_sormlq.c core_blas/core_sormqr.c core_blas/core_spamm.c
  core_blas/core_spemv.c core_blas/core_spotrf.c core_blas/core_ssygst.c
  core_blas/core_ssymm.c core_blas/core_ssyr2k.c core_blas/core_ssyrk.c
  core_blas/core_ssyssq.c core_blas/core_stradd.c core_blas/core_strmm.c
  core_blas/core_strssq.c core_blas/core_strtri.c core_blas/core_stslqt.c
  core_blas/core_stsmlq.c core_blas/core_stsmqr.c core_blas/core_stsqrt.c
  core_blas/core_sttlqt.c core_blas/core_sttmlq.c core_blas/core_sttmqr.c
  core_blas/core_sttqrt.c
)

set(TEST_SOURCE
  test/test.h test/test.c
  test/test_dzamax.c test/test_damax.c test/test_scamax.c test/test_samax.c
  test/test_zcposv.c test/test_dsposv.c test/test_zgbsv.c test/test_dgbsv.c
  test/test_cgbsv.c test/test_sgbsv.c test/test_zgbtrf.c test/test_dgbtrf.c
  test/test_cgbtrf.c test/test_sgbtrf.c test/test_zgeadd.c test/test_dgeadd.c
  test/test_cgeadd.c test/test_sgeadd.c test/test_zgelqf.c test/test_dgelqf.c
  test/test_cgelqf.c test/test_sgelqf.c test/test_zgelqs.c test/test_dgelqs.c
  test/test_cgelqs.c test/test_sgelqs.c test/test_zgels.c test/test_dgels.c
  test/test_cgels.c test/test_sgels.c test/test_zgemm.c test/test_dgemm.c
  test/test_cgemm.c test/test_sgemm.c test/test_zgeqrf.c test/test_dgeqrf.c
  test/test_cgeqrf.c test/test_sgeqrf.c test/test_zgeqrs.c test/test_dgeqrs.c
  test/test_cgeqrs.c test/test_sgeqrs.c
  test/test_zcgbsv.c test/test_dsgbsv.c test/test_zgetri.c test/test_dgetri.c
  test/test_cgetri.c test/test_sgetri.c test/test_zgetri_aux.c
  test/test_dgetri_aux.c test/test_cgetri_aux.c test/test_sgetri_aux.c
  test/test_zhemm.c test/test_chemm.c test/test_zher2k.c test/test_cher2k.c
  test/test_zherk.c test/test_cherk.c test/test_zhetrf.c test/test_dsytrf.c
  test/test_chetrf.c test/test_ssytrf.c test/test_zhesv.c test/test_dsysv.c
  test/test_chesv.c test/test_ssysv.c test/test_zlacpy.c test/test_dlacpy.c
  test/test_clacpy.c test/test_slacpy.c test/test_zlag2c test/test_clag2z.c
  test/test_dlag2s.c test/test_slag2d.c test/test_zlange.c test/test_dlange.c
  test/test_clange.c test/test_slange.c test/test_zlanhe.c test/test_clanhe.c
  test/test_zlansy.c test/test_dlansy.c test/test_clansy.c test/test_slansy.c
  test/test_zlantr.c test/test_dlantr.c test/test_clantr.c test/test_slantr.c
  test/test_zlascl.c test/test_dlascl.c test/test_clascl.c test/test_slascl.c
  test/test_zlaset.c test/test_dlaset.c test/test_claset.c test/test_slaset.c
  test/test_zgeswp.c test/test_dgeswp.c test/test_cgeswp.c test/test_sgeswp.c
  test/test_zlauum.c test/test_dlauum.c test/test_clauum.c test/test_slauum.c
  test/test_zpbsv.c test/test_dpbsv.c test/test_cpbsv.c test/test_spbsv.c
  test/test_zpbtrf.c test/test_dpbtrf.c test/test_cpbtrf.c test/test_spbtrf.c
  test/test_zlangb.c test/test_dlangb.c test/test_clangb.c test/test_slangb.c
  test/test_zposv.c test/test_dposv.c test/test_cposv.c test/test_sposv.c
  test/test_zpoinv.c test/test_dpoinv.c test/test_cpoinv.c test/test_spoinv.c
  test/test_zpotrf.c test/test_dpotrf.c test/test_cpotrf.c test/test_spotrf.c
  test/test_zpotri.c test/test_dpotri.c test/test_cpotri.c test/test_spotri.c
  test/test_zpotrs.c test/test_dpotrs.c test/test_cpotrs.c test/test_spotrs.c
  test/test_zsymm.c test/test_dsymm.c test/test_csymm.c test/test_ssymm.c
  test/test_zsyr2k.c test/test_dsyr2k.c test/test_csyr2k.c test/test_ssyr2k.c
  test/test_zsyrk.c test/test_dsyrk.c test/test_csyrk.c test/test_ssyrk.c
  test/test_ztradd.c test/test_dtradd.c test/test_ctradd.c test/test_stradd.c
  test/test_ztrmm.c test/test_dtrmm.c test/test_ctrmm.c test/test_strmm.c
  test/test_ztrsm.c test/test_dtrsm.c test/test_ctrsm.c test/test_strsm.c
  test/test_ztrtri.c test/test_dtrtri.c test/test_ctrtri.c test/test_strtri.c
  test/test_zunmlq.c test/test_dormlq.c test/test_cunmlq.c test/test_sormlq.c
  test/test_zunmqr.c test/test_dormqr.c test/test_cunmqr.c test/test_sormqr.c
)

add_library(plasma OBJECT
  include/plasma.h
  ${COMPUTE_SOURCE}
  control/tree.c control/tuning.c control/workspace.c control/version.c
)

# CMake knows about "plasma" library at this point so inform CMake where the headers are
target_include_directories(plasma PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

add_library(plasma_core_blas OBJECT include/plasma_core_blas.h
  ${CORE_SOURCE}
  control/barrier.c control/async.c
)

target_include_directories(plasma_core_blas PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

add_executable(plasmatest
  include/plasma.h
  ${TEST_SOURCE}
)

find_library(MATH_LIBRARY m)
if( MATH_LIBRARY )
  # OpenBLAS needs to link C math library (usually -lm) but MKL doesn't
  set(PLASMA_LIBRARIES ${PLASMA_LINALG_LIBRARIES} ${LUA_LIBRARIES} ${MATH_LIBRARY})
else( MATH_LIBRARY )
  set(PLASMA_LIBRARIES ${PLASMA_LINALG_LIBRARIES} ${LUA_LIBRARIES})
endif( MATH_LIBRARY )

if(USE_TAGGING)
  set(PLASMA_LIBRARIES ${PLASMA_LIBRARIES} omptracing)
endif()

target_link_libraries( plasmatest plasma plasma_core_blas ${PLASMA_LIBRARIES} )
target_link_libraries( plasma plasma_core_blas ${PLASMA_LIBRARIES} )
target_link_libraries( plasma_core_blas ${PLASMA_LIBRARIES} )

target_include_directories(plasmatest PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

configure_file( include/plasma_config.hin ${CMAKE_CURRENT_SOURCE_DIR}/include/plasma_config.h @ONLY NEWLINE_STYLE LF )

install(TARGETS plasma plasma_core_blas LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
file( GLOB plasma_headers include/plasma*.h)
install(FILES ${plasma_headers} DESTINATION include)
install(TARGETS plasmatest RUNTIME DESTINATION bin)
