noinst_PROGRAMS=fwdpy11_cpp_tests

TESTS=$(noinst_PROGRAMS)

fwdpy11_cpp_tests_SOURCES=fwdpy11_cpp_tests.cc \
						  test_discrete_demography_stable_sorting_of_events.cc \
						  test_lowlevel_demographic_events.cc \
						  test_lowlevel_discrete_demography_bad_models.cc \
						  test_lowlevel_discrete_demography_objects.cc \
						  discrete_demography_roundtrips.cc \
						  discrete_demography_util.cc \
						  test_MutationDominance.cc 

AM_CPPFLAGS=-I../fwdpy11/headers -I../fwdpy11/headers/fwdpp -I../fwdpy11/src/evolve_population
AM_CXXFLAGS=-W -Wall --coverage -DBOOST_TEST_DYN_LINK

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

coverage-local: check coverage.sh
	sh coverage.sh

clean-local:
	find . -name '*.gcov' | xargs rm -f
	find . -name '*.gcda' | xargs rm -f
	find . -name '*.gcno' | xargs rm -f
	find ../fwdpy11/src -name '*.gcov' | xargs rm -f
	find ../fwdpy11/src -name '*.gcda' | xargs rm -f
	find ../fwdpy11/src -name '*.gcno' | xargs rm -f
	rm -f *.info
	rm -rf fwdpy11_coverage

