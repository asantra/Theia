
include_directories(Tracker)
include_directories(${THEIAPATH}/Common/Particles/source/Particles)
include_directories(${THEIAPATH}/Common/Magnets/source/Magnets)
include_directories(${THEIAPATH}/Common/Logger/source/Logger)

file(GLOB SRC
    "src/*.cxx"
)

root_generate_dictionary(Tracker_Dict TrkDetector.h Material.h Cluster.h TrkLayer.h TrkProbe.h DrawProbe.h TrackParameter.h TTreeStream.h LINKDEF Tracker/Tracker_LinkDef.h)
add_library(Tracker SHARED ${SRC} Tracker_Dict.cxx)
#### adding executables
add_executable(TrackReconstruction Root/Reconstruction.cxx Tracker_Dict.cxx)
add_executable(TrackDigitization Root/Digitization.cxx Tracker_Dict.cxx)



target_link_libraries(Tracker PUBLIC ${ROOT_LIBRARIES} ROOT::Geom PUBLIC Magnets Logger)
target_include_directories(Tracker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Tracker)

### for executables
target_link_libraries(TrackReconstruction PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)
target_include_directories(TrackReconstruction PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)

target_link_libraries(TrackDigitization PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)
target_include_directories(TrackDigitization PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)



install(TARGETS Tracker DESTINATION ${CMAKE_SOURCE_DIR}/../install/lib)
### for executables
install(TARGETS TrackReconstruction DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)
install(TARGETS TrackDigitization DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)


















