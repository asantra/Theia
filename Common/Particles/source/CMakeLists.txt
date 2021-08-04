

include_directories(Particles)

file(GLOB SRC
    "src/*.cxx"
)

root_generate_dictionary(Particles_Dict Particles.h LUXEUtility.h LINKDEF Particles/Particles_LinkDef.h)
add_library(Particles SHARED ${SRC} Particles_Dict.cxx)

target_link_libraries(Particles PUBLIC ${ROOT_LIBRARIES})

target_include_directories (Particles PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


install(TARGETS Particles DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
