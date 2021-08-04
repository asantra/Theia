
include_directories(Tracker)
include_directories(${THEIAPATH}/Common/Particles/source/Particles)
include_directories(${THEIAPATH}/Common/Magnets/source/Magnets)
include_directories(${THEIAPATH}/Common/Logger/source/Logger)

file(GLOB SRC
    "src/*.cxx"
)

root_generate_dictionary(Tracker_Dict TrkDetector.h Tracker.h Material.h Cluster.h TrkLayer.h TrkProbe.h DrawProbe.h TrackParameter.h TTreeStream.h LINKDEF Tracker/Tracker_LinkDef.h)
add_library(Tracker SHARED ${SRC} Tracker_Dict.cxx)

target_link_libraries(Tracker PUBLIC ${ROOT_LIBRARIES} ROOT::Geom PUBLIC Magnets Logger)

target_include_directories(Tracker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Tracker)


install(TARGETS Tracker DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)


