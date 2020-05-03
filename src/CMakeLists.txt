
add_library(physica STATIC
    Linear_Algebra.cpp
    Natural_Units.cpp
    Numerics.cpp
    Statistics.cpp
    Utilities.cpp
)

target_compile_options(physica PUBLIC -Wall -pedantic)

target_link_libraries(physica
	PUBLIC coverage_config 
	)

target_include_directories(physica  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS physica DESTINATION ${LIB_DIR})
