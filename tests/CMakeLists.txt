# 1. Linear algebra
add_executable(test_Linear_Algebra test_Linear_Algebra.cpp)
target_link_libraries(test_Linear_Algebra 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Linear_Algebra PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Linear_Algebra PUBLIC -Wall -pedantic)
install(TARGETS test_Linear_Algebra DESTINATION ${TESTS_DIR})
add_test(NAME Test_Linear_Algebra COMMAND test_Linear_Algebra
	WORKING_DIRECTORY ${TESTS_DIR})

# 2. Natural Units
add_executable(test_Natural_Units test_Natural_Units.cpp)
target_link_libraries(test_Natural_Units 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Natural_Units PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Natural_Units PUBLIC -Wall -pedantic)
install(TARGETS test_Natural_Units DESTINATION ${TESTS_DIR})
add_test(NAME Test_Natural_Units COMMAND test_Natural_Units
	WORKING_DIRECTORY ${TESTS_DIR})

# 3. Numerics
add_executable(test_Numerics test_Numerics.cpp)
target_link_libraries(test_Numerics 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Numerics PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Numerics PUBLIC -Wall -pedantic)
install(TARGETS test_Numerics DESTINATION ${TESTS_DIR})
add_test(NAME Test_Numerics COMMAND test_Numerics
	WORKING_DIRECTORY ${TESTS_DIR})

# 4. Statistics
add_executable(test_Statistics test_Statistics.cpp)
target_link_libraries(test_Statistics 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Statistics PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Statistics PUBLIC -Wall -pedantic)
install(TARGETS test_Statistics DESTINATION ${TESTS_DIR})
add_test(NAME Test_Statistics COMMAND test_Statistics
	WORKING_DIRECTORY ${TESTS_DIR})

# 5. Utilities
add_executable(test_Utilities test_Utilities.cpp)
target_link_libraries(test_Utilities 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Utilities PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Utilities PUBLIC -Wall -pedantic)
install(TARGETS test_Utilities DESTINATION ${TESTS_DIR})
add_test(NAME Test_Utilities COMMAND test_Utilities
	WORKING_DIRECTORY ${TESTS_DIR})