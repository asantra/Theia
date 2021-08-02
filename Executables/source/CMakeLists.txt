cmake_minimum_required(VERSION 3.10)

project(Executables LANGUAGES CXX)
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)

### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()

add_executable(TrackReconstruction Executables/Root/Reconstruction.cxx)

target_link_libraries(TrackReconstruction PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)

target_include_directories(TrackReconstruction PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)


# ### installation
install(TARGETS TrackReconstruction DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)
