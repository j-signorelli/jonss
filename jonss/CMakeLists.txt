# Add the relative path (from main directory) to the list of file names
# (Inspired largely from convert_filenames_to_full_paths in MFEM)
# (Copied directly from https://github.com/j-signorelli/jabber/blob/main/jabber/CMakeLists.txt)
function(append_relative_paths NAMES)
   unset(tmp_names)
   foreach(name ${${NAMES}})
      set(rel_path "")
      cmake_path(RELATIVE_PATH CMAKE_CURRENT_SOURCE_DIR
                  BASE_DIRECTORY ${CMAKE_SOURCE_DIR} 
                  OUTPUT_VARIABLE rel_path)
      list(APPEND tmp_names ${rel_path}/${name})
   endforeach()
   set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()

set(JONSS_SOURCES
      integrators.cpp)

set(JONSS_HEADERS
      options.hpp
      constants.hpp
      state.hpp
      primitives.hpp
      fluxes.hpp
      integrators.hpp
      jonss.hpp)

append_relative_paths(JONSS_SOURCES)
append_relative_paths(JONSS_HEADERS)

set(JONSS_SOURCES ${JONSS_SOURCES} PARENT_SCOPE)
set(JONSS_HEADERS ${JONSS_HEADERS} PARENT_SCOPE)