cmake_minimum_required(VERSION 3.14)
project(bv_tree)

#### USER DEFINED ##############################################################

#### BASIC SETTINGS ############################################################

include_directories(.)

set (CMAKE_CXX_FLAGS "-std=c++2a  -msse4.2 -Wall -Wextra -O3 -g") #with assertions


#google test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#### TARGETS ###################################################################

add_executable(main src/main.cpp)



#### TESTS ###################################################################
enable_testing()

add_executable(hello_test test/hello_test.cpp)
target_link_libraries(hello_test gtest_main)

add_executable(test_naive_bv test/test_naive_bv.cpp)
target_link_libraries(test_naive_bv gtest_main)

add_executable(test_naive_tree_bp test/test_naive_tree_bp.cpp)
target_link_libraries(test_naive_tree_bp gtest_main)

add_executable(test_static_small_bv test/test_static_small_bv.cpp)
target_link_libraries(test_static_small_bv gtest_main)

add_executable(test_dynamic_bv test/test_dynamic_bv.cpp)
target_link_libraries(test_dynamic_bv gtest_main)

add_executable(test_dynamic_tree_bp test/test_dynamic_tree_bp.cpp)
target_link_libraries(test_dynamic_tree_bp gtest_main)

include(GoogleTest)
gtest_discover_tests(hello_test)
