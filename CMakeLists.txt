cmake_minimum_required(VERSION 3.10)

project(Theia LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)


# if (UNIX)
#   # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c -O3 -fPIC -Wall -Wextra -std=c++11 -m64")
#   # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -fPIC -m64 -std=gnu++17")
# endif()

add_compile_options(-std=c++17)

find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


# Add your dependencies in the order you want them to build.
add_subdirectory(Common/ConfigReader/source)
add_subdirectory(Common/Logger/source)
add_subdirectory(Common/Magnets/source)
add_subdirectory(Common/Particles/source)
add_subdirectory(Common/Plotting/source)
add_subdirectory(Tracker/source)
add_subdirectory(Executables/source)