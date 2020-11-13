noinst_PROGRAMS=fwdpy11_cpp_tests

TESTS=$(noinst_PROGRAMS)

fwdpy11_cpp_tests_SOURCES=fwdpy11_cpp_tests.cc \
						  test_discrete_demography_stable_sorting_of_events.cc

AM_CPPFLAGS=-I../fwdpy11/headers -I../fwdpy11/headers/fwdpp 
AM_CXXFLAGS=-W -Wall --coverage -DBOOST_TEST_DYN_LINK

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

coverage-local: check coverage.sh
	sh coverage.sh

clean-local:
	find . -name '*.gcov' | xargs rm -f
	find . -name '*.gcda' | xargs rm -f
	find . -name '*.gcno' | xargs rm -f
	rm -f *.info
	rm -rf fwdpy11_coverage

