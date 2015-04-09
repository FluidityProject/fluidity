Point(1) = {0, 0, 2.7, 1.0};
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{2};
}
Extrude {{-1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{-1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{5};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{4, 3, 1, 2};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{8, 5, 12, 16};
}
Extrude {{0, -1, 0}, {0, 0, 0}, Pi/4} {
  Line{4, 3, 1, 2};
}
Extrude {{0, -1, 0}, {0, 0, 0}, Pi/4} {
  Line{33, 36, 40, 44};
}
Transfinite Line {47} = 1 Using Progression 1;
Transfinite Line {33} = 1 Using Progression 1;
Transfinite Line {4}  = 1 Using Progression 1;
Transfinite Line {5}  = 1 Using Progression 1;
Transfinite Line {23} = 1 Using Progression 1;
Transfinite Line {50} = 1 Using Progression 1;
Transfinite Line {36} = 1 Using Progression 1;
Transfinite Line {3}  = 1 Using Progression 1;
Transfinite Line {8}  = 1 Using Progression 1;
Transfinite Line {19} = 1 Using Progression 1;
Transfinite Line {54} = 1 Using Progression 1;
Transfinite Line {40} = 1 Using Progression 1;
Transfinite Line {1}  = 1 Using Progression 1;
Transfinite Line {12} = 1 Using Progression 1;
Transfinite Line {26} = 1 Using Progression 1;
Transfinite Line {58} = 1 Using Progression 1;
Transfinite Line {44} = 1 Using Progression 1;
Transfinite Line {2}  = 1 Using Progression 1;
Transfinite Line {16} = 1 Using Progression 1;
Transfinite Line {30} = 1 Using Progression 1;
Transfinite Line {48} = 1 Using Progression 1;
Transfinite Line {34} = 1 Using Progression 1;
Transfinite Line {6}  = 1 Using Progression 1;
Transfinite Line {21} = 1 Using Progression 1;
Transfinite Line {51} = 1 Using Progression 1;
Transfinite Line {37} = 1 Using Progression 1;
Transfinite Line {9}  = 1 Using Progression 1;
Transfinite Line {20} = 1 Using Progression 1;
Transfinite Line {56} = 1 Using Progression 1;
Transfinite Line {42} = 1 Using Progression 1;
Transfinite Line {14} = 1 Using Progression 1;
Transfinite Line {28} = 1 Using Progression 1;
// top
Transfinite Surface {49} Alternated;
Transfinite Surface {35} Alternated;
Transfinite Surface {7} Alternated;
Transfinite Surface {25} Alternated;
// upper middle
Transfinite Surface {53} Alternated;
Transfinite Surface {39} Alternated;
Transfinite Surface {11} Alternated;
Transfinite Surface {22} Alternated;
// lower middle
Transfinite Surface {57} Alternated;
Transfinite Surface {43} Alternated;
Transfinite Surface {15} Alternated;
Transfinite Surface {29} Alternated;
// bottom
Transfinite Surface {60} Alternated;
Transfinite Surface {46} Alternated;
Transfinite Surface {18} Alternated;
Transfinite Surface {32} Alternated;
// flat edges
Physical Line(17) = {47};
Physical Line(18) = {23};
Physical Line(19) = {50};
Physical Line(20) = {19};
Physical Line(21) = {54};
Physical Line(22) = {26};
Physical Line(23) = {58};
Physical Line(24) = {30};
// top
Physical Surface(33) = {49};
Physical Surface(34) = {35};
Physical Surface(35) = {7};
Physical Surface(36) = {25};
// upper middle
Physical Surface(25) = {53};
Physical Surface(26) = {39};
Physical Surface(27) = {11};
Physical Surface(28) = {22};
// lower middle
Physical Surface(29) = {57};
Physical Surface(30) = {43};
Physical Surface(31) = {15};
Physical Surface(32) = {29};
// bottom
Physical Surface(37) = {60};
Physical Surface(38) = {46};
Physical Surface(39) = {18};
Physical Surface(40) = {32};
