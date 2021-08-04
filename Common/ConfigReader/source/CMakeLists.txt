
include_directories(ConfigReader)

file(GLOB SRC
    "src/*.cxx"
)

### generate root dictionary, essential for building shared library
root_generate_dictionary(ConfigReader_Dict ConfigReader.h LINKDEF ConfigReader/ConfigReader_LinkDef.h)
add_library(ConfigReader SHARED ${SRC} ConfigReader_Dict.cxx)


target_link_libraries(ConfigReader PUBLIC ${ROOT_LIBRARIES})

target_include_directories (ConfigReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ConfigReader DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
