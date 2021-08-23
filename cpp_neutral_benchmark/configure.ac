AC_INIT([cpp_neutral_benchmark], 0.0)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)
AC_CONFIG_FILES(Makefile)

dnl check for things that are required to compile all examples and/or use library at all
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(boost/program_options.hpp, ,[AC_MSG_ERROR([echo "boost/program_options.hpp not found."])])

dnl check for C run-time libraries
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])

dnl check for C++ runtime libraries
AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
	  AC_CHECK_LIB([boost_program_options],[main],
				   HAVEBOOSTPORUNTIME=1,
				   [AC_MSG_ERROR([boost_program_options runtime library does not seem to be here.])])
      AC_REQUIRE_CPP
AC_LANG_RESTORE
AC_OUTPUT
