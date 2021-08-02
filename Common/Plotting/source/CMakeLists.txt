cmake_minimum_required(VERSION 3.10)
# Create a library called "Tracker" which includes the source files in src/ and inc/.
# The extension is already found. Any number of sources could be listed here.

project(Hists LANGUAGES CXX)

# Finding the ROOT package
find_package(ROOT CONFIG REQUIRED)

# arka+noam
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


include_directories(Plotting)
#include_directories(src)

file(GLOB SRC
    "src/*.cxx"
)

add_library(Plotting SHARED ${SRC})
target_link_libraries(Plotting PUBLIC ROOT::RIO ROOT::Net ROOT::Physics ROOT::Core ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::Rint ROOT::MathMore ROOT::Graf3d)

# Make sure the compiler can find include files for our Tracker library
# when other libraries or executables link to Tracker
target_include_directories (Plotting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ROOT::RIO ROOT::Net ROOT::Physics ROOT::Core ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::Rint ROOT::MathMore ROOT::Graf3d)

### installation
install(TARGETS Plotting DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
