add_subdirectory(Catch2)
add_executable(unit-tests
                test_primitives.cpp
                test_fluxes.cpp)

target_link_libraries(unit-tests PRIVATE Catch2::Catch2WithMain jonss)

include(CTest)
include(Catch)
catch_discover_tests(unit-tests)