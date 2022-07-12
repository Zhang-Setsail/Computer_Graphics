cmake_minimum_required(VERSION 3.10)
project(HW2)

find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)

add_executable(HW2 main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
target_link_libraries(HW2 ${OpenCV_LIBRARIES})
