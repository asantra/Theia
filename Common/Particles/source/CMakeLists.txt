cmake_minimum_required(VERSION 3.10)

project(ParticleClass LANGUAGES CXX)
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)

### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()

include_directories(Particles)

file(GLOB SRC
    "src/*.cxx"
)

ROOT_GENERATE_DICTIONARY(Particles_Dict Particles.h LUXEUtility.h LINKDEF Particles/Particles_LinkDef.h)
add_library(Particles SHARED ${SRC} Particles_Dict.cxx)

target_link_libraries(Particles PUBLIC ${ROOT_LIBRARIES})

target_include_directories (Particles PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories (Particles PUBLIC ${THEIAPATH}/Common/Particles/source)


install(TARGETS Particles DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
