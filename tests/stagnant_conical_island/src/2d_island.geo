dx = 3.75;
Point(1) = {0.0,-4.45,0, dx};
Point(2) = {0.0,-0.85,0, dx};
Line(1) = {2,1};
Point(3) = {0.81240384046359604,0,0, dx};
Extrude {0,-0.25,0} {
  Point{3};
}
Extrude {3.6,0,0} {
  Line{2};
}
Rotate {{0,0,1}, {0,0,0}, Pi/14} {
  Duplicata { Line{1}; }
}
Rotate {{0,0,1}, {0,0,0}, Pi/14} {
  Duplicata { Line{7}; }
}
Rotate {{0,0,1}, {0,0,0}, Pi/14} {
  Duplicata { Line{8}; }
}
Rotate {{0,0,1}, {0,0,0}, Pi/14} {
  Duplicata { Line{9}; }
}
Rotate {{0,0,1}, {0,0,0}, Pi/13} {
  Duplicata { Line{10}; }
}
Extrude {3.4537158056565742,-1.0157987663695807,0} {
  Point{4};
}

Line(13) = {1,8};
Line(14) = {8,10};
Line(15) = {10,12};
Line(16) = {12,14};
Line(17) = {14,16};
Line(18) = {16,17};
Line(19) = {17,6};
Line(20) = {17,15};
Line(21) = {15,4};
Line(22) = {15,13};
Line(23) = {13,11};
Line(24) = {11,9};
Line(25) = {9,7};
Line(26) = {7,2};
Extrude {-0.81240384046359604,0.0,0} {
  Line{2};
}
Point(20) = {0.0,-0.60104076400856532,0, dx};

Line(31) = {19,20};
Line(32) = {20,4};
Line(33) = {20,2};
Line Loop(34) = {1,13,-7,26};
Plane Surface(35) = {34};
Line Loop(36) = {7,14,-8,25};
Plane Surface(37) = {36};
Line Loop(38) = {15,-9,24,8};
Plane Surface(39) = {38};
Line Loop(40) = {9,16,-10,23};
Plane Surface(41) = {40};
Line Loop(42) = {10,17,-11,22};
Plane Surface(43) = {42};
Line Loop(44) = {11,18,20};
Plane Surface(45) = {44};
Line Loop(46) = {20,21,12};
Plane Surface(47) = {46};
Line Loop(48) = {12,19,-5};
Plane Surface(49) = {48};
Line Loop(50) = {25,26,-33,32,-21,22,23,24};
Plane Surface(51) = {50};
Line Loop(52) = {32,29,31};
Plane Surface(53) = {52};
Point(21) = {0.0,-7.5,0, dx};
Point(22) = {12.5,-7.5,0, dx};
Point(23) = {12.5,0,0, dx};
Line(54) = {21,22};
Line(55) = {22,23};
Line(56) = {23,5};
Line(57) = {1,21};
Line Loop(66) = {54,55,56,3,-19,-18,-17,-16,-15,-14,-13,57};
Plane Surface(67) = {66};
Physical Line(68) = {54};
Physical Line(69) = {55};
Physical Line(70) = {56, 4, 28};
Physical Line(71) = {27, 31, 33, 1, 57};
Physical Surface(72) = {67};
Physical Surface(73) = {35};
Physical Surface(74) = {37};
Physical Surface(75) = {39};
Physical Surface(76) = {41};
Physical Surface(77) = {43};
Physical Surface(78) = {45};
Physical Surface(79) = {47};
Physical Surface(80) = {49};
Physical Surface(81) = {6};
Physical Surface(82) = {51};
Physical Surface(83) = {53};
Physical Surface(84) = {30};

Field[99] = MathEval;
Field[99].F = "3.75";

// 81MainRamp
Field[81] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[81].F = "1.0-(1.0-0.01)*exp(-((0.33333333333333337*x - 0.0*y + 1.529198719845468 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[181] = Restrict;
Field[181].IField = 81;
Field[181].FacesList = {6};
Field[181].EdgesList = {4, 5};

// 80TopTriangle
Field[80] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[80].F = "1.0-(1.0-0.01)*exp(-((0.33333333333333326*x - 0.04800300977795665*y + 1.517197967400979 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[180] = Restrict;
Field[180].IField = 80;
Field[180].FacesList = {49};
Field[180].EdgesList = {5, 12};

// 79MiddleTriangle
Field[79] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[79].F = "1.0-(1.0-0.01)*exp(-((0.31234419682042946*x - 0.11936607392182969*y + 1.5164088564761085 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[179] = Restrict;
Field[179].IField = 79;
Field[179].FacesList = {47};
Field[179].EdgesList = {12, 20};

// 78LowerTriangle
Field[78] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[78].F = "1.0-(1.0-0.01)*exp(-((0.31277924212070413*x - 0.11769616749659653*y + 1.5166666666666671 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[178] = Restrict;
Field[178].IField = 78;
Field[178].FacesList = {45};
Field[178].EdgesList = {11, 20};

// 77IslandSideFifthLeftThin
Field[77] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[77].F = "1.0-(1.0-0.01)*exp(-((0.28584561690952376*x - 0.17618609090395509*y + 1.5166666666666668 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[177] = Restrict;
Field[177].IField = 77;
Field[177].FacesList = {43};
Field[177].EdgesList = {10, 11};

// 76IslandSideFourthLeft
Field[76] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[76].F = "1.0-(1.0-0.01)*exp(-((0.2371936844982894*x - 0.23719368449828956*y + 1.5166666666666668 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[176] = Restrict;
Field[176].IField = 76;
Field[176].FacesList = {41};
Field[176].EdgesList = {9, 10};

// 75IslandSideThirdLeft
Field[75] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[75].F = "1.0-(1.0-0.01)*exp(-((0.17846618340753292*x - 0.28402730381373031*y + 1.5166666666666668 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[175] = Restrict;
Field[175].IField = 75;
Field[175].FacesList = {39};
Field[175].EdgesList = {8, 9};

// 74IslandSideSecondLeft
Field[74] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[74].F = "1.0-(1.0-0.01)*exp(-((0.11078964267083921*x - 0.31661860812121567*y + 1.5166666666666668 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[174] = Restrict;
Field[174].IField = 74;
Field[174].FacesList = {37};
Field[174].EdgesList = {7, 8};

// 73IslandSideLeft
Field[73] = MathEval;
//         maxele - (maxele-minele)*exp(-((mx*x + my*y + c - h)^2)/(2*sigma^2))
//         h = 3.0 - ((4.45 - R0)/3.0)      
Field[73].F = "1.0-(1.0-0.01)*exp(-((0.037557646633370328*x - 0.33333333333333331*y + 1.5166666666666666 - 2.2066666666666666)^2)/(2*0.05^2))";
Field[173] = Restrict;
Field[173].IField = 73;
Field[173].FacesList = {35};
Field[173].EdgesList = {1, 7};

Field[100] = Min;
Field[100].FieldsList = {99, 173, 174, 175, 176, 177, 178, 179, 180, 181};

Background Field = 100;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

