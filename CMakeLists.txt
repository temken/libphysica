cmake_minimum_required(VERSION 3.12)

project(template_cpp_cmake)

# Use C++-11
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

set(BIN_DIR     ${PROJECT_SOURCE_DIR}/bin)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR     ${PROJECT_SOURCE_DIR}/lib)
set(SRC_DIR     ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR   ${PROJECT_SOURCE_DIR}/tests)
set(EXTERNAL_DIR  ${PROJECT_SOURCE_DIR}/external)

include_directories( ${INCLUDE_DIR} )

add_subdirectory( ${SRC_DIR} )

include(FetchContent)
# Googletest
set(GTEST_DIR   ${EXTERNAL_DIR}/googletest)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
  SOURCE_DIR     "${GTEST_DIR}/src"
  BINARY_DIR     "${GTEST_DIR}/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)

# if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)
# endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")


option (BUILD_TESTING "Build the testing tree." ON)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
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