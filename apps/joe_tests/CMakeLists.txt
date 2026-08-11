
# Require jabber
find_package(jabber REQUIRED)

add_executable(jonss_wave_test jonss_wave_test.cpp utils.cpp)
target_link_libraries(jonss_wave_test PRIVATE mfem jabber jonss)

install(TARGETS jonss_wave_test 
   DESTINATION ${CMAKE_INSTALL_BINDIR})