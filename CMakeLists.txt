cmake_minimum_required(VERSION 3.17)

project(libphysica VERSION "0.1.4")

# Use C++-14
set(CMAKE_CXX_STANDARD 14)

# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

set(GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external)

# External projects
include(ExternalProject)
include(FetchContent)

# 1. Boost
find_package(Boost 1.65 REQUIRED)

# 2. libconfig++
option(DOWNLOAD_LIBCONFIG "Download the libconfig++ library without checking for installed versions" OFF)

if(NOT DOWNLOAD_LIBCONFIG)
  find_path(LIBCONFIG_INCLUDE_DIRs libconfig.h++
    /usr/local/include
    /usr/include
    ~/libs/include
    ~/include)
  find_library(LIBCONFIGPP_LIBRARY_PATH config++
    /usr/local/lib
    /usr/lib
    ~/libs/lib
    ~/lib)

  if(LIBCONFIGPP_LIBRARY_PATH AND LIBCONFIG_INCLUDE_DIRs)
    message(STATUS "Found libconfig++: ${LIBCONFIGPP_LIBRARY_PATH} and headers in ${LIBCONFIG_INCLUDE_DIRs}.")
  else()
    message(WARNING "libconfig++ could not be found. The library will be downloaded.")
    set(DOWNLOAD_LIBCONFIG ON)
  endif()
endif()

if(DOWNLOAD_LIBCONFIG)
  message(STATUS "libconfig++: Download version 1.7.2 from http://www.hyperrealm.com/libconfig/ and install it in /external/.")
  ExternalProject_Add(libconfig_download
    PREFIX ${EXTERNAL_DIR}/libconfig
    URL https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ${EXTERNAL_DIR}/libconfig/src/libconfig_download/configure
    --prefix=${EXTERNAL_DIR}/libconfig
    --enable-shared
    --disable-examples
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
  )
  set(LIBCONFIG_INCLUDE_DIRs ${EXTERNAL_DIR}/libconfig/include)
  set(LIBCONFIGPP_LIBRARY_PATH ${EXTERNAL_DIR}/libconfig/lib/libconfig++${CMAKE_SHARED_LIBRARY_SUFFIX})
  message(STATUS "libconfig++ will be installed in ${LIBCONFIGPP_LIBRARY_PATH} and headers in ${LIBCONFIG_INCLUDE_DIRs}.")
endif()

add_library(libconfig++ SHARED IMPORTED)
set_target_properties(libconfig++ PROPERTIES IMPORTED_LOCATION ${LIBCONFIGPP_LIBRARY_PATH})

# target_include_directories(libconfig++ INTERFACE ${LIBCONFIG_INCLUDE_DIRs})

# 3. Googletest
set(GTEST_DIR ${EXTERNAL_DIR}/googletest)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0
  SOURCE_DIR "${GTEST_DIR}/src"
  BINARY_DIR "${GTEST_DIR}/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  TEST_COMMAND ""
)

# Version header
execute_process( # Git commit hash
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process( # Git commit hash
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process( # Git commit hash
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CMAKE_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process( # Git commit hash
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CMAKE_GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
  ${INCLUDE_DIR}/version.hpp.in
  ${GENERATED_DIR}/version.hpp
)

# Source and include directories
include_directories(${INCLUDE_DIR} ${GENERATED_DIR})
add_subdirectory(${SRC_DIR})

# Code Coverage Configuration
if(NOT TARGET coverage_config)
  add_library(coverage_config INTERFACE)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)

if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0 # no optimization
    -g # generate debug info
    --coverage # sets all required flags
  )

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

option(BUILD_TESTING "Build the testing tree." ON)

# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if(BUILD_TESTING AND(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()

  FetchContent_GetProperties(googletest)

  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(
      ${googletest_SOURCE_DIR}
      ${googletest_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif()

  add_subdirectory(${TESTS_DIR})
endif()