
add_executable(tests
	test_Linear_Algebra.cpp)
target_link_libraries(tests 
	PRIVATE
		physica
		gtest_main	#contains the main function
)
# coverage_config is contained in module_library.

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR})
