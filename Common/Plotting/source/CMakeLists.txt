cmake_minimum_required(VERSION 3.10)

project(Hists LANGUAGES CXX)

# Finding the ROOT package
find_package(ROOT CONFIG REQUIRED)


include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


include_directories(Plotting)

file(GLOB SRC
    "src/*.cxx"
)

add_library(Plotting SHARED ${SRC})
target_link_libraries(Plotting PUBLIC ROOT::RIO ROOT::Net ROOT::Physics ROOT::Core ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::Rint ROOT::MathMore ROOT::Graf3d)


target_include_directories (Plotting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ROOT::RIO ROOT::Net ROOT::Physics ROOT::Core ROOT::Graf ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::Rint ROOT::MathMore ROOT::Graf3d)

### installation
install(TARGETS Plotting DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
