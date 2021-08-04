
include_directories(Magnets)


file(GLOB SRC
    "src/*.cxx"
)

### generate root dictionary, essential for building shared library
root_generate_dictionary(Magnets_Dict Magnets.h Utils.h LINKDEF Magnets/Magnets_LinkDef.h)

add_library(Magnets SHARED ${SRC} Magnets_Dict.cxx)

target_include_directories(Magnets PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Magnets PUBLIC ${ROOT_LIBRARIES} ROOT::Geom ${ROOT_EXE_LINKER_FLAGS})

install(TARGETS Magnets DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
