cmake_minimum_required(VERSION 3.15)

# Test executable
add_executable(run_tests test_main.cpp)

target_link_libraries(run_tests PRIVATE option_pricer_core)

target_include_directories(run_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Add test
add_test(NAME OptionPricerTests COMMAND run_tests)
