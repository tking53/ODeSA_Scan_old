
add_executable(scan Scan.cpp)
target_include_directories(scan PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(scan ${ROOT_LIBRARIES})
install(TARGETS scan DESTINATION ${CMAKE_INSTALL_PREFIX})
