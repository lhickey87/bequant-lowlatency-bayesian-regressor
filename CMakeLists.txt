cmake_minimum_required(VERSION 3.20)
project(MyProject VERSION 0.1 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

add_library(mylib STATIC
    src/OLS.cpp
    src/Matrix.cpp
)

target_include_directories(mylib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

add_executable(OLS
    main.cpp
)

target_link_libraries(OLS PRIVATE
    mylib
    OpenMP::OpenMP_CXX
)
