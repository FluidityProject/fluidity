include_directories(. JudyCommon Judy1 JudyL)
add_definitions(-DJU_64BIT)
add_subdirectory(JudyCommon)
add_subdirectory(Judy1)
add_subdirectory(JudyL)
add_subdirectory(JudySL)
add_subdirectory(JudyHS)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Judy1/Judy1Tables.c
	COMMAND Judy1TablesGen
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Judy1)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/JudyL/JudyLTables.c
	COMMAND JudyLTablesGen
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/JudyL)

add_library(judy STATIC ${JUDY_SRCS} ${JUDY_SRCS_JUDY1}
                 ${JUDY_SRCS_JUDYL} ${JUDY_SRCS_JUDYSL} ${JUDY_SRCS_JUDYH}
		 ${CMAKE_CURRENT_BINARY_DIR}/Judy1/Judy1Tables.c
		 ${CMAKE_CURRENT_BINARY_DIR}/JudyL/JudyLTables.c)

set_source_files_properties(${JUDY_SRCS_JUDY1} PROPERTIES
	COMPILE_DEFINITIONS JUDY1)
set_source_files_properties(${JUDY_SRCS_JUDYL} PROPERTIES
	COMPILE_DEFINITIONS JUDYL)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/Judy1/Judy1Tables.c
	PROPERTIES COMPILE_DEFINITIONS JUDY1)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/JudyL/JudyLTables.c
	PROPERTIES COMPILE_DEFINITIONS JUDYL)

set_source_files_properties(Judy1/j__udy1Test.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDY1 -DJUDYGETINLINE")
set_source_files_properties(Judy1/Judy1Next.c Judy1/Judy1NextEmpty.c PROPERTIES
	COMPILE_FLAGS  "${CMAKE_C_FLAGS} -DJUDY1 -DJUDYNEXT")
set_source_files_properties(Judy1/Judy1Prev.c Judy1/Judy1PrevEmpty.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDY1 -DJUDYPREV")
set_source_files_properties(Judy1/Judy1ByCount.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDY1 -DNOSMARTJBB -DNOSMARTJBU -DNOSMARTJLB")

set_source_files_properties(JudyL/j__udyLTest.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDYL -DJUDYGETINLINE")
set_source_files_properties(JudyL/JudyLNext.c JudyL/JudyLNextEmpty.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDYL -DJUDYNEXT")
set_source_files_properties(JudyL/JudyLPrev.c JudyL/JudyLPrevEmpty.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDYL -DJUDYPREV")
set_source_files_properties(JudyL/JudyLByCount.c PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDYL -DNOSMARTJBB -DNOSMARTJBU -DNOSMARTJLB")

install(TARGETS judy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
