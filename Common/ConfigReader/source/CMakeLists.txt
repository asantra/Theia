cmake_minimum_required(VERSION 3.10)

project(ConfigReaderClass LANGUAGES CXX)
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


### get the path of THEIA
if (DEFINED ENV{THEIAPATH})
  set(THEIAPATH "$ENV{THEIAPATH}" CACHE INTERNAL "Copied from environment variable")
endif()


include_directories(ConfigReader)

file(GLOB SRC
    "src/*.cxx"
)

ROOT_GENERATE_DICTIONARY(ConfigReader_Dict ConfigReader.h LINKDEF ConfigReader/ConfigReader_LinkDef.h)
add_library(ConfigReader SHARED ${SRC} ConfigReader_Dict.cxx)


target_link_libraries(ConfigReader PUBLIC ${ROOT_LIBRARIES})

target_include_directories (ConfigReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ConfigReader DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
