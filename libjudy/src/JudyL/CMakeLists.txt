add_executable(JudyLTablesGen EXCLUDE_FROM_ALL JudyLTablesGen.c)
set_target_properties(JudyLTablesGen PROPERTIES
	COMPILE_DEFINITIONS JUDYL)

set(JUDY_SRCS_JUDYL JudyL/JudyLByCount.c JudyL/JudyLCount.c
		    JudyL/JudyLDecascade.c JudyL/JudyLFirst.c
    		    JudyL/j__udyLGet.c JudyL/JudyLInsArray.c
		    JudyL/JudyLInsertBranch.c JudyL/JudyLMemActive.c
		    JudyL/JudyLNext.c JudyL/JudyLPrev.c JudyL/JudyLCascade.c
		    JudyL/JudyLCreateBranch.c JudyL/JudyLDel.c
    		    JudyL/JudyLFreeArray.c JudyL/JudyLGet.c JudyL/JudyLIns.c
		    JudyL/JudyLMallocIF.c JudyL/JudyLMemUsed.c
		    JudyL/JudyLNextEmpty.c JudyL/JudyLPrevEmpty.c
       PARENT_SCOPE)
