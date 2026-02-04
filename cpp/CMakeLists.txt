cmake_minimum_required(VERSION 3.15)
project(OptionPricerCore VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Create library
add_library(option_pricer_core
    src/models/BlackScholesModel.cpp
    src/products/VanillaOption.cpp
    src/products/AsianOption.cpp
    src/products/AmericanOption.cpp
    src/products/BarrierOption.cpp
    src/pricing/AnalyticPricer.cpp
    src/pricing/MonteCarloPricer.cpp
    src/pricing/MCSettings.cpp
    src/pricing/FDSettings.cpp
    src/pricing/FDPricer.cpp
    src/greeks/BSGreeks.cpp
    src/utils/MathUtils.cpp
    src/utils/ImpliedVol.cpp
    src/utils/FDUtils.cpp
)

target_include_directories(option_pricer_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Link math library on Unix systems
if(UNIX AND NOT APPLE)
    target_link_libraries(option_pricer_core PUBLIC m)
endif()

# Enable testing
enable_testing()
add_subdirectory(tests)

# Installation rules
install(TARGETS option_pricer_core
    EXPORT OptionPricerCoreTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include)
