add_executable(fldecomp fldgmsh.cpp fldmain.cpp fldtriangle.cpp partition.cpp)
target_link_libraries(fldecomp fluidity)
