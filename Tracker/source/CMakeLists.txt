cmake_minimum_required(VERSION 3.10)

project(TrackerClass LANGUAGES CXX)

find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()


include_directories(Tracker)
include_directories(${THEIAPATH}/Common/Particles/source/Particles)
include_directories(${THEIAPATH}/Common/Magnets/source/Magnets)
include_directories(${THEIAPATH}/Common/Logger/source/Logger)

file(GLOB SRC
    "src/*.cxx"
)

ROOT_GENERATE_DICTIONARY(Tracker_Dict Detector.h Tracker.h Material.h Cluster.h Layer.h Probe.h DrawProbe.h TrackParameter.h TTreeStream.h LINKDEF Tracker/Tracker_LinkDef.h)
add_library(Tracker SHARED ${SRC} Tracker_Dict.cxx)

target_link_libraries(Tracker PUBLIC ${ROOT_LIBRARIES} ROOT::Geom PUBLIC Magnets Logger)

target_include_directories(Tracker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Tracker)

target_include_directories(Tracker PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets)
target_include_directories(Tracker PUBLIC ${THEIAPATH}/Common/Logger/source/Logger)
target_include_directories(Tracker PUBLIC ${THEIAPATH}/Common/Particles/source/Particles)


install(TARGETS Tracker DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)


