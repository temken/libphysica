
add_library(libphysica STATIC
    Linear_Algebra.cpp
    Natural_Units.cpp
    Numerics.cpp
    Statistics.cpp
    Utilities.cpp
)

set_target_properties(libphysica PROPERTIES PREFIX "")

target_compile_options(libphysica PUBLIC -Wall -pedantic)

target_link_libraries(libphysica
	PUBLIC coverage_config 
	)

target_include_directories(libphysica
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		${INCLUDE_DIR})

install(TARGETS libphysica DESTINATION ${LIB_DIR})
