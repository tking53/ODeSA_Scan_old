project(Scan_E20003)

cmake_minimum_required(VERSION 3.0)

#------------------------------------------------------------------
# Assign install prefix to `install` if not defined.
#
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	message(STATUS "Install Prefix not specified.")
	get_filename_component(INSTALL_DIR ${PROJECT_SOURCE_DIR}/install REALPATH)
	set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH
			"Install Prefix" FORCE)
endif ()
message(STATUS "Installing to ${CMAKE_INSTALL_PREFIX}")

#----------------------------------------------------------------------------
# Define the default build type
#
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Build type, options are: None Debug Release
			RelWithDebInfo MinSizeRel."
		FORCE)
	message(STATUS "Build type not defined, using default.")
ENDIF(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#----------------------------------------------------------------------------
# Treat warnings as errors for testing.
#
option(WARNINGS_AS_ERRORS "Treat all warnings as errors." OFF)
mark_as_advanced(WARNINGS_AS_ERRORS)

if (WARNINGS_AS_ERRORS)
	message(STATUS "Warnings are being treated as errors.")
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -pedantic -Wall -Wextra")
	else ()
		message(WARNING "Unknown compiler to set warnings as errors.")
	endif ()
endif(WARNINGS_AS_ERRORS)


find_package(ROOT REQUIRED COMPONENTS Spectrum)

add_subdirectory(src)


