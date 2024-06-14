project(Cubic)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package (Eigen3 REQUIRED NO_MODULE)

add_executable(cubic main.cpp)

include_directories(./include)
 
target_link_libraries(cubic PUBLIC Eigen3::Eigen)

target_link_libraries(cubic PUBLIC ncurses)


