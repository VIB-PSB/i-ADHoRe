find_package(GoogleTest 1.0.0)

if (GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})
    add_executable(test PackingTest.cpp test.cpp
        indexToXYTest.cpp ParToolBoxTest.cpp
        IntroduceGapsTest.cpp AlignTest.cpp ../src/higherLevel.cpp
        ../src/threadPool.cpp
        ../src/parallel.cpp ../src/levelTwo.cpp ../src/BaseCluster.cpp
        ../src/Cluster.cpp ../src/DataSet.cpp ../src/GHM.cpp
        ../src/GHMProfile.cpp ../src/Gene.cpp ../src/GeneFamily.cpp
        ../src/GenePairs.cpp ../src/GeneList.cpp ../src/AlignmentDrawer
        ../src/ListElement.cpp ../src/Multiplicon.cpp ../src/SvgWriter.cpp
        ../src/Profile.cpp ../src/Settings.cpp ../src/hpmath.cpp
        ../src/util.cpp ../src/alignComp.cpp ../src/SynthenicCloud.cpp)
    target_link_libraries(test ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})
    target_link_libraries(test ${MPI_LIBRARIES})
    target_link_libraries(test ${CMAKE_THREAD_LIBS_INIT})
    target_link_libraries(test bmp alignment png)
endif(GTEST_FOUND)

# Setup the variables we use below.
# set(PARSER_TESTS
#   test/ws_parser_test.cpp
#   test/label_node_test.cpp
#   test/symbol_node_test.cpp
#   test/immediate_number_node_test.cpp
#   test/immediate_node_test.cpp
#   test/internal_name_node_test.cpp
#   test/dereferenced_name_node_test.cpp
#   test/scoped_name_node_test.cpp
#   test/indexed_name_node_test.cpp
#   test/compound_name_node_test.cpp
#   test/name_node_test.cpp
#   test/expression_node_test.cpp
#   test/operation_node_test.cpp
#   test/function_call_node_test.cpp
#   test/procedure_call_node_test.cpp
#   test/mutation_node_test.cpp
#   test/statement_node_test.cpp
#   test/if_statement_node_test.cpp
#   test/block_entry_node_test.cpp
#   test/block_node_test.cpp)

# if(GTEST_FOUND)
#   include_directories(${GTEST_INCLUDE_DIRS})
#   foreach(TEST ${PARSER_TESTS})
#     string(REGEX REPLACE "\\.(cc|cxx|cpp)" "" TEST_NAME ${TEST})
#     string(REPLACE "test/" "" TEST_NAME ${TEST_NAME})
#     set(TEST_SOURCE ${TEST})
#     add_executable(${TEST_NAME} ${TEST_SOURCE})
#     set_target_properties(${TEST_NAME} PROPERTIES
#       RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")
#     target_link_libraries(${TEST_NAME} inc_parser
#       ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})
#   endforeach(TEST)
# endif(GTEST_FOUND)
