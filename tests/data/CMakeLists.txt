foreach(mesh IN LISTS unittest_meshes_3d)
    add_custom_target(${mesh}
    	COMMAND gmsh -3 ${mesh}.geo
	WORKING_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})
endforeach()

foreach(mesh IN LISTS unittest_meshes_2d)
    add_custom_target( ${mesh}
    	COMMAND gmsh -2 ${mesh}.geo
	WORKING_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})
endforeach()
