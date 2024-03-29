cmake_minimum_required(VERSION 3.2)
project(traycer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT LINUX AND NOT APPLE)
    message(FATAL_ERROR "Target system is not supported.")
endif()

add_compile_options(-Wall -Wextra -Wpedantic -DGLM_FORCE_RADIANS)

add_subdirectory(vendor/glm)

add_library(cli cli.cpp)

add_library(scene scene.cpp)
target_link_libraries(scene PUBLIC glm)
target_include_directories(scene PUBLIC vendor)

add_executable(traycer main.cpp)
target_link_libraries(traycer PRIVATE glm cli scene)
target_include_directories(traycer PRIVATE vendor)

if(LINUX)
    target_link_libraries(traycer PRIVATE -lGLEW -lGL -lglut)
elseif(APPLE)
    target_link_libraries(traycer PRIVATE "-framework OpenGL" "-framework GLUT")
    target_compile_options(traycer PRIVATE -Wno-deprecated-declarations)
endif()