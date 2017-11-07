PROJECT(PolynomialConvexCones)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

#eigen
set(CMAKE_MODULE_PATH $ENV{Eigen3_DIR}/cmake)
set(EIGEN3_INCLUDE_DIR $ENV{Eigen3_DIR})
FIND_PACKAGE(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#qhull
set(QHULL_PATH $ENV{QHULL_ROOT})
include_directories(${QHULL_PATH}/include)
link_directories(${QHULL_PATH}/lib)

ADD_EXECUTABLE(
	sample
	sample.cpp
	PCCCalculator.cpp
	PCCCalculator.h
	Combination.cpp
	Combination.h
)

TARGET_LINK_LIBRARIES(sample qhullstatic.lib qhullstatic_r.lib)	
