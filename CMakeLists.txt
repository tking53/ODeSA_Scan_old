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
find_package(ROOT REQUIRED COMPONENTS Spectrum)

add_subdirectory(src)


