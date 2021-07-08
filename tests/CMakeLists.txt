# 1. Integration
add_executable(test_Integration test_Integration.cpp)
target_link_libraries(test_Integration 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Integration PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Integration PUBLIC -Wall -pedantic)
install(TARGETS test_Integration DESTINATION ${TESTS_DIR})
add_test(NAME Test_Integration COMMAND test_Integration
	WORKING_DIRECTORY ${TESTS_DIR})

# 2. Linear algebra
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

# 3. Natural Units
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

# 4. Numerics
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

# 5. Special Functions
add_executable(test_Special_Functions test_Special_Functions.cpp)
target_link_libraries(test_Special_Functions 
	PRIVATE
		libphysica
		gtest_main	#contains the main function
)
target_include_directories(test_Special_Functions PRIVATE ${GENERATED_DIR} )
target_compile_options(test_Special_Functions PUBLIC -Wall -pedantic)
install(TARGETS test_Special_Functions DESTINATION ${TESTS_DIR})
add_test(NAME Test_Special_Functions COMMAND test_Special_Functions
	WORKING_DIRECTORY ${TESTS_DIR})

# 6. Statistics
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

# 7. Utilities
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