cmake_minimum_required(VERSION 3.10)

project(LoggerClass LANGUAGES CXX)
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)

### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()

include_directories(Logger)

file(GLOB SRC
    "src/*.cxx"
)

ROOT_GENERATE_DICTIONARY(Logger_Dict LUXELog.h LINKDEF Logger/LUXELog_LinkDef.h)
add_library(Logger SHARED ${SRC} Logger_Dict.cxx)

target_link_libraries(Logger PUBLIC ${ROOT_LIBRARIES})


target_include_directories(Logger PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


### installation
install(TARGETS Logger DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)

