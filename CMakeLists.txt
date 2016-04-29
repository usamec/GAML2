cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11 -O2 -pg")

enable_testing()
find_package(GTest REQUIRED)
find_package (Threads REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_library(node node.cc)

add_executable(node_test node_test.cc)
target_link_libraries(node_test node ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(NodeTest node_test)

add_library(graph graph.cc)
target_link_libraries(graph node)

add_executable(graph_test graph_test.cc)
target_link_libraries(graph_test graph ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(GraphTest graph_test)

add_library(path path.cc)
target_link_libraries(path node)
add_executable(path_test path_test.cc)
target_link_libraries(path_test path ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} graph)
add_test(PathTest path_test)

add_library(read_set read_set.cc)
add_executable(read_set_test read_set_test.cc)
target_link_libraries(read_set_test read_set ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(ReadSetTest read_set_test)

add_executable(util_test util_test.cc)
target_link_libraries(util_test ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(UtilTest util_test)

add_library(path_aligner path_aligner.cc)
target_link_libraries(path_aligner path read_set graph)

add_library(moves moves.cc)
target_link_libraries(moves path)

add_executable(moves_test moves_test.cc)
target_link_libraries(moves_test ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} moves)
add_test(MovesTest moves_test)

add_library(read_probability_calculator read_probability_calculator.cc)
target_link_libraries(read_probability_calculator path_aligner)

add_executable(read_probability_calculator_test read_probability_calculator_test.cc)
target_link_libraries(read_probability_calculator_test read_probability_calculator ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(ReadProbabilityCalculatorTest read_probability_calculator_test)

add_executable(get_subgraph get_subgraph.cc)
target_link_libraries(get_subgraph graph)

add_executable(gaml gaml_main.cc)
target_link_libraries(gaml graph path read_probability_calculator moves)

