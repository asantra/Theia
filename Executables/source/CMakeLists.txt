

add_executable(TrackReconstruction Executables/Root/Reconstruction.cxx)

target_link_libraries(TrackReconstruction PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)

target_include_directories(TrackReconstruction PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)


# ### installation
install(TARGETS TrackReconstruction DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)
