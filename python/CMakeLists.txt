cmake_minimum_required(VERSION 3.15)
project(OptionPricerPython)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    # Try to find it via Python
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE pybind11_RESULT
    )
    if(pybind11_RESULT EQUAL 0)
        find_package(pybind11 CONFIG REQUIRED PATHS ${pybind11_DIR})
    else()
        message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
    endif()
endif()

# Add the C++ core library
add_subdirectory(../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp_build)

# Create Python module
pybind11_add_module(option_pricer_cpp bindings.cpp)

target_link_libraries(option_pricer_cpp PRIVATE option_pricer_core)

# Static link MinGW runtime on Windows to avoid DLL dependency issues
if(MINGW)
    target_link_options(option_pricer_cpp PRIVATE 
        -static-libgcc 
        -static-libstdc++ 
        -static
        -lpthread
    )
endif()

# Installation
install(TARGETS option_pricer_cpp LIBRARY DESTINATION .)
