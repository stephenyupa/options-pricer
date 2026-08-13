cmake_minimum_required(VERSION 3.15)
project(options_pricer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_TESTS "Build unit tests" ON)

add_compile_options(-Wall -Wextra)

add_library(options_pricer STATIC
    src/normal_dist.cpp
    src/black_scholes.cpp
    src/monte_carlo.cpp
    src/binomial_tree.cpp
)
target_include_directories(options_pricer PUBLIC include)

add_executable(validation_harness app/main.cpp)
target_link_libraries(validation_harness PRIVATE options_pricer)

if(BUILD_TESTS)
    enable_testing()

    add_executable(test_black_scholes tests/test_black_scholes.cpp)
    target_link_libraries(test_black_scholes PRIVATE options_pricer)
    add_test(NAME test_black_scholes COMMAND test_black_scholes)

    add_executable(test_binomial_tree tests/test_binomial_tree.cpp)
    target_link_libraries(test_binomial_tree PRIVATE options_pricer)
    add_test(NAME test_binomial_tree COMMAND test_binomial_tree)

    add_executable(test_monte_carlo tests/test_monte_carlo.cpp)
    target_link_libraries(test_monte_carlo PRIVATE options_pricer)
    add_test(NAME test_monte_carlo COMMAND test_monte_carlo)

    add_executable(test_parity tests/test_parity.cpp)
    target_link_libraries(test_parity PRIVATE options_pricer)
    add_test(NAME test_parity COMMAND test_parity)
endif()
