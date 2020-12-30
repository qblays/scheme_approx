#cmake_minimum_required(VERSION 3.10)
#project(matplot_test)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O2 ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -O2 ")

#find_package(Matplot++ REQUIRED)

#add_executable(matplot_test plot_example.cpp)
#target_link_libraries(matplot_test PUBLIC Matplot++::matplot)
cmake_minimum_required(VERSION 3.10)
project(prakt)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O2 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Ofast -march=native -g")

find_package(Matplot++ REQUIRED)

add_executable(prakt main.cpp task_1.cpp task_2.cpp task_3.cpp task_4.cpp common.cpp)
target_link_libraries(prakt PUBLIC Matplot++::matplot boost_program_options)
install(TARGETS prakt
	RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}")
