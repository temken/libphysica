file(GLOB FILES "*.cpp")
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(libphysica STATIC
        ${FILES} )

set_target_properties(libphysica PROPERTIES PREFIX "")

target_compile_options(libphysica PUBLIC -Wall -pedantic)

# Find and include libconfig
find_path(LIBCONFIG_INCLUDE_DIRs libconfig.h++
    /usr/local/include
    /usr/include
    ~/libs/include
    ~/include )
find_library(LIBCONFIGPP_LIBRARY config++
    /usr/local/lib
    /usr/lib
    ~/libs/lib
    ~/lib )

target_link_libraries(libphysica
    PRIVATE ${Boost_LIBRARIES}
	PUBLIC
        coverage_config
        ${LIBCONFIGPP_LIBRARY} 
    )

target_include_directories(libphysica
	PRIVATE
		${GENERATED_DIR}
        ${Boost_INCLUDE_DIR}
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		${INCLUDE_DIR}
        ${LIBCONFIG_INCLUDE_DIRs} 
        )

install(TARGETS libphysica DESTINATION ${LIB_DIR})
