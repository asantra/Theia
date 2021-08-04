
include_directories(Logger)

file(GLOB SRC
    "src/*.cxx"
)

### generate root dictionary, essential for building shared library
root_generate_dictionary(Logger_Dict LUXELog.h LINKDEF Logger/LUXELog_LinkDef.h)
add_library(Logger SHARED ${SRC} Logger_Dict.cxx)

target_link_libraries(Logger PUBLIC ${ROOT_LIBRARIES})


target_include_directories(Logger PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


### installation
install(TARGETS Logger DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)

