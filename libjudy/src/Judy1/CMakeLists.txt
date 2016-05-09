include_directories(.)
add_executable(Judy1TablesGen EXCLUDE_FROM_ALL Judy1TablesGen.c)
set_target_properties(Judy1TablesGen PROPERTIES
	COMPILE_FLAGS "${CMAKE_C_FLAGS} -DJUDY1")

set(JUDY_SRCS_JUDY1 Judy1/Judy1Test.c
	     Judy1/j__udy1Test.c Judy1/Judy1Set.c Judy1/Judy1SetArray.c
	     Judy1/Judy1Unset.c Judy1/Judy1First.c Judy1/Judy1Next.c 
    	     Judy1/Judy1Prev.c Judy1/Judy1NextEmpty.c Judy1/Judy1PrevEmpty.c
	     Judy1/Judy1Count.c Judy1/Judy1ByCount.c Judy1/Judy1FreeArray.c
	     Judy1/Judy1MemUsed.c Judy1/Judy1MemActive.c Judy1/Judy1Cascade.c
	     Judy1/Judy1Decascade.c Judy1/Judy1CreateBranch.c
	     Judy1/Judy1InsertBranch.c Judy1/Judy1MallocIF.c
	     PARENT_SCOPE)
