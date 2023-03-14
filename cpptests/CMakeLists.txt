SET(CPPTEST_SOURCES
    fwdpy11_cpp_tests.cc
    test_MutationDominance.cc
    test_gss_classes.cc
    test_effect_sizes_of_zero.cc
    test_Sregion_from_mvnorm.cc
    test_forward_demes_graph.cc
    forward_demes_graph_fixtures.cc
    test_evolvets.cc
    test_core_genetic_map_regions.cc
)

add_executable(fwdpy11_cpp_tests ${CPPTEST_SOURCES})
add_dependencies(fwdpy11_cpp_tests fwdpy11core header)
target_link_directories(fwdpy11_cpp_tests PRIVATE ${CMAKE_SOURCE_DIR}/fwdpy11)
target_link_libraries(fwdpy11_cpp_tests PRIVATE boost_unit_test_framework fwdpy11core GSL::gsl GSL::gslcblas)
if (NOT APPLE)
    target_link_options(fwdpy11_cpp_tests BEFORE PUBLIC LINKER:--no-as-needed -ldl)
endif()
set_target_properties(fwdpy11_cpp_tests PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_TEST_DYN_LINK")

add_test(NAME fwdpy11_cpp_tests
    COMMAND fwdpy11_cpp_tests)
