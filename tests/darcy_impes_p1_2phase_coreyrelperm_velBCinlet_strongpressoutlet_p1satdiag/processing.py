import darcy_impes_tools

numericalFilenameStem = "darcy_impes_p1_2phase_coreyrelperm_velBCinlet_strongpressoutlet_p1satdiag"
modelNameList = ["relpermupwind", "modrelpermupwind_satfesweby"]
gridNameListPerDimension = [["A", "B", "C", "D"], 
                            ["A", "B", "C"],
                            ["A", "B"]]
fieldNameList = ["Phase2::Saturation"]
analyticFilenameStem = "reference_solution/analytic_BL_CoreyPerm"

testHelper = darcy_impes_tools.TestHelper(numericalFilenameStem, modelNameList, gridNameListPerDimension, analyticFilenameStem, fieldNameList)
testHelper.processFolder()
testHelper.generateReports()

