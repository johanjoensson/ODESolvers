cmake_minimum_required(VERSION 3.14)

project(
    ODE
	VERSION 0.1
    DESCRIPTION "C++ template library for ODE solvers"
	LANGUAGES CXX
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

add_library(ode-lib INTERFACE)

target_include_directories(ode-lib INTERFACE 
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	)

if(BUILD_PYTHON)
	find_package(pybind11)
	if(pybind11_FOUND)
		add_subdirectory(python)
	else()
		message(STATUS "pybind11 not found, not building python bindings.")
	endif()
endif()

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(CMAKE_CXX_EXTENSIONS OFF)
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	find_package(Doxygen)
	if(Doxygen_FOUND)
		add_subdirectory(docs)
	else()
		message(STATUS "Doxygen not found, not building docs")
	endif()

	set(WFLAGS -Werror -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor
		-Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic
		-Wconversion -Wsign-conversion -Wmisleading-indentation
		-Wnull-dereference -Wdouble-promotion -Wformat=2 -Weffc++)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set(WFLAGS ${WFLAGS} -Wduplicated-cond -Wduplicated-branches
			-Wlogical-op -Wuseless-cast)
	endif()
endif() # Main project only part ends here

# IDEs should put the headers in a nice place
source_group(
	TREE "${PROJECT_SOURCE_DIR}/include"
	PREFIX "Header Files"
	FILES ${HEADER_LIST}
)

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
AND BUILD_TESTING)
    enable_testing()
	add_subdirectory(tests)
endif()

install(TARGETS ode-lib 
	EXPORT ode-targets
	INCLUDES DESTINATION include/ode)
install(EXPORT ode-targets
	NAMESPACE ODE::
	DESTINATION share/ode
	)
install(DIRECTORY include/ DESTINATION include/ode 
	FILES_MATCHING PATTERN "*.h")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ode/ode-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT ode-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ode/ode-targets.cmake"
  NAMESPACE ODE::
)
configure_file(cmake/ode-config.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/ode/ode-config.cmake"
  COPYONLY
)

set(ConfigPackageLocation share/cmake/ode)
install(EXPORT ode-targets
  FILE
    ode-targets.cmake
  NAMESPACE
  ODE::
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    cmake/ode-config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/ode/ode-config-version.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
