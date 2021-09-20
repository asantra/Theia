
#include_directories(Executables)

#root_generate_dictionary(Executables_Dict Root/Executables.h LINKDEF Root/Executables_LinkDef.h)

add_executable(TrackReconstruction Executables/Root/Reconstruction.cxx Executables/Root/Executables.h)
add_executable(TrackDigitization Executables/Root/Digitization.cxx Executables/Root/Executables.h)

target_link_libraries(TrackReconstruction PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)
target_include_directories(TrackReconstruction PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)


target_link_libraries(TrackDigitization PUBLIC ${ROOT_LIBRARIES} PUBLIC Tracker Logger Magnets Particles)
target_include_directories(TrackDigitization PUBLIC ${THEIAPATH}/Common/Magnets/source/Magnets ${THEIAPATH}/Common/Particles/source/Particles ${THEIAPATH}/Tracker/source/Tracker ${THEIAPATH}/Common/Logger/source/Logger ${THEIAPATH}/Common/Plotting/source/Plotting)

# ### installation
install(TARGETS TrackReconstruction DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)
install(TARGETS TrackDigitization DESTINATION ${CMAKE_SOURCE_DIR}/../install/bin)
