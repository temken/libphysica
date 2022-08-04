file(GLOB FILES "*.cpp")
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(libphysica STATIC
    ${FILES})

if(DOWNLOAD_LIBCONFIG)
    add_dependencies(libphysica libconfig_download) # ensures that libconfig gets build before libphysica.
endif()

set_target_properties(libphysica PROPERTIES PREFIX "")

target_compile_options(libphysica PUBLIC -Wall -pedantic)

target_link_libraries(libphysica
    PRIVATE ${Boost_LIBRARIES}
    PUBLIC
    coverage_config
    libconfig++
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
