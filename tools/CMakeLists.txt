file(GLOB local_scripts RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py *.sh)

foreach(script IN LISTS local_scripts)
    get_filename_component(script_name ${script} NAME_WE)
    list(APPEND scripts ${script_name})
    add_custom_target( ${script_name} ALL
          "${CMAKE_COMMAND}" "-E" copy "${CMAKE_CURRENT_SOURCE_DIR}/${script}" ${CMAKE_BINARY_DIR}/bin/${script_name})
endforeach()
