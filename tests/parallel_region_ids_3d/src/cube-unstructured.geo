layers=2;
Point(1) = {0, 0, 0, 0.25};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{2}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{1}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{4}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{2}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{6}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{3}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{8}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{5}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{7}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{4}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{6}; Layers{layers};
}
Extrude {0, 0, -1} {
  Point{1}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{2}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{3}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{8}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{6}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{4}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{5}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{7}; Layers{2*layers};
}
Extrude {0, 0, -1} {
  Point{9}; Layers{2*layers};
}
Extrude {0.5, 0, 0} {
  Point{10}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{11}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{15}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{14}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{16}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{17}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{10}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{15}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{11}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{14}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{12}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{13}; Layers{layers};
}
Line Loop(34) = {23, 32, -25, -30};
Plane Surface(35) = {34};
Line Loop(36) = {22, 30, -24, -28};
Plane Surface(37) = {36};
Line Loop(38) = {24, 31, -26, -29};
Plane Surface(39) = {38};
Line Loop(40) = {25, 33, -27, -31};
Plane Surface(41) = {40};
Line Loop(42) = {25, -16, -12, 17};
Plane Surface(43) = {42};
Line Loop(44) = {17, -24, -18, 11};
Plane Surface(45) = {44};
Line Loop(46) = {30, -17, -5, 14};
Plane Surface(47) = {46};
Line Loop(48) = {6, 20, -31, -17};
Plane Surface(49) = {48};
Line Loop(50) = {26, -20, -9, 19};
Plane Surface(51) = {50};
Line Loop(52) = {27, -21, -10, 20};
Plane Surface(53) = {52};
Line Loop(54) = {33, -21, -8, 16};
Plane Surface(55) = {54};
Line Loop(56) = {32, -16, -7, 15};
Plane Surface(57) = {56};
Line Loop(58) = {13, 28, -18, -3};
Plane Surface(59) = {58};
Line Loop(60) = {18, 29, -19, -4};
Plane Surface(61) = {60};
Line Loop(62) = {22, -14, -1, 13};
Plane Surface(63) = {62};
Line Loop(64) = {23, -15, -2, 14};
Plane Surface(65) = {64};
Line Loop(66) = {2, 7, -12, -5};
Plane Surface(67) = {66};
Line Loop(68) = {12, 8, -10, -6};
Plane Surface(69) = {68};
Line Loop(70) = {11, 6, -9, -4};
Plane Surface(71) = {70};
Line Loop(72) = {11, -5, -1, 3};
Plane Surface(73) = {72};
Surface Loop(74) = {37, 63, 73, 59, 45, 47};
Volume(75) = {74};
Surface Loop(76) = {35, 65, 57, 67, 47, 43};
Volume(77) = {76};
Surface Loop(78) = {55, 41, 53, 69, 49, 43};
Volume(79) = {78};
Surface Loop(80) = {51, 39, 61, 71, 45, 49};
Volume(81) = {80};
Physical Surface(82) = {65, 63};
Physical Surface(83) = {59, 61};
Physical Surface(84) = {51, 53};
Physical Surface(85) = {55, 57};
Physical Surface(86) = {69, 71, 73, 67};
Physical Surface(87) = {41, 35, 39, 37};
Physical Volume(28) = {81};
Physical Volume(27) = {79};
Physical Volume(25) = {75};
Physical Volume(26) = {77};
