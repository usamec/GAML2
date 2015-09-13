cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")

enable_testing()
find_package(GTest REQUIRED)
find_package (Threads REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_library(node node.cc)

add_executable(node_test node_test.cc)
target_link_libraries(node_test node ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(NodeTest node_test)

add_library(graph graph.cc)

add_executable(graph_test graph_test.cc)
target_link_libraries(graph_test graph ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_test(GraphTest graph_test)