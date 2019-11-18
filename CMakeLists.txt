cmake_minimum_required(VERSION 3.12)
project(cpp_finance_week7_code)

set(CMAKE_CXX_STANDARD 14)

add_executable(cpp_finance_week7_code main.cpp MonteCarlo.cpp MonteCarlo.hpp Random.cpp Random.hpp Utility.cpp Utility.hpp)