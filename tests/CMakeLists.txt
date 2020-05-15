
add_executable(tests
	test_Linear_Algebra.cpp
	test_Natural_Units.cpp
	test_Numerics.cpp
	test_Statistics.cpp
	test_Utilities.cpp)

target_link_libraries(tests 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
# coverage_config is contained in module_library.

target_include_directories(tests PRIVATE ${GENERATED_DIR} )

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR})
