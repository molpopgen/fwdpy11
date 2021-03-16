AC_INIT([fwdpy11_cpp_tests], 0.0)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)
AC_CONFIG_FILES(Makefile)

dnl check for things that are required to compile all examples and/or use library at all
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(boost/test/unit_test.hpp, ,[AC_MSG_ERROR([echo "boost/test/unit_test.hpp not found. Unit tests will not
											be compiled."])])

dnl check for C run-time libraries
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])

dnl check for C++ runtime libraries
AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_REQUIRE_CPP
      AC_CHECK_LIB(boost_unit_test_framework,main,HAVEBOOSTTESTRUNTIME=1,[AC_MSG_ERROR([boost_unit_test_framework runtime library does not
                                                            seem to be here.  Unit tests will not be compiled])])
AC_LANG_RESTORE
AC_OUTPUT
