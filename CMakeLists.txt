cmake_minimum_required(VERSION 3.15)

project(jonss LANGUAGES CXX CUDA VERSION 2.0.0)

# Get MFEM
find_package(MFEM 4.9 REQUIRED CONFIG) 

#add_library(jonss)
#target_compile_features(jonss PUBLIC cxx_std_20) # Require C++20
#add_subdirectory(jonss)
# target_sources(jonss
#    PRIVATE 
#       ${JONSS_SOURCES}
#    PUBLIC
#       FILE_SET public_headers TYPE HEADERS
#       FILES
#          ${JONSS_HEADERS})


add_executable(run_jonss run_jonss.cpp)
target_link_libraries(run_jonss PRIVATE mfem)

# Install jonss
include(GNUInstallDirs)

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# From: https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake
install(TARGETS run_jonss
   EXPORT jonssTargets
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT jonssTargets
   FILE "jonss-targets.cmake"
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jonss)
   
include(CMakePackageConfigHelpers)
configure_package_config_file("cmake/jonss-config.cmake.in"
   "jonss-config.cmake"
   INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jonss)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/jonss-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES 
   ${CMAKE_CURRENT_BINARY_DIR}/jonss-config.cmake
   ${CMAKE_CURRENT_BINARY_DIR}/jonss-config-version.cmake
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/jonss)
