cmake_minimum_required(VERSION 3.10)

project(Theia LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)


### this is not necessary, so commented out
# add_compile_options(-std=c++17)

find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)

### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()


# Add your dependencies in the order you want them to build.
add_subdirectory(Common/ConfigReader/source)
add_subdirectory(Common/Logger/source)
add_subdirectory(Common/Magnets/source)
add_subdirectory(Common/Particles/source)
add_subdirectory(Common/Plotting/source)
add_subdirectory(Tracker/source)

