cmake_minimum_required(VERSION 3.10)

project(MagnetClass LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()

include_directories(Magnets)


file(GLOB SRC
    "src/*.cxx"
)

ROOT_GENERATE_DICTIONARY(Magnets_Dict Magnets.h Utils.h LINKDEF Magnets/Magnets_LinkDef.h)
add_library(Magnets SHARED ${SRC} Magnets_Dict.cxx)

target_include_directories(Magnets PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Magnets PUBLIC ${ROOT_LIBRARIES} ROOT::Geom ${ROOT_EXE_LINKER_FLAGS})

install(TARGETS Magnets DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
