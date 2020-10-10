cmake_minimum_required(VERSION 3.15)
project(NewUnderwaterSFM)

set(CMAKE_CXX_STANDARD 11)

add_executable(NewUnderwaterSFM main.cpp UnderwaterSFM.cpp UnderwaterSFM.h)

#set(CMAKE_CXX_FLAGS "-std=c++11 -O0")   #original sign is -O2, change it to -O0 can fix the "optimized out" bug
set(PCL_DIR "/usr/share/pcl-1.8/")

find_package(OpenCV 3.1 REQUIRED)
find_package(Ceres REQUIRED)
find_package(Sophus REQUIRED)
find_package(PCL REQUIRED COMPONENT common io)
include_directories(${PCL_INCLUDE_DIRS})
include_directories( ${Sophus_INCLUDE_DIRS} )
include_directories("/usr/include/eigen3")
add_definitions(${PCL_DEFINITIONS})

target_link_libraries(NewUnderwaterSFM ${OpenCV_LIBS} ${CERES_LIBRARIES} ${Sophus_LIBRARIES} ${PCL_LIBRARIES})


