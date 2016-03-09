//Earth's Radius
radius = 6.37101e+06;
//typical mesh element size
lc=radius/Pi;

//Point at the4 centre of the planet
Point(1) = {0.0, 0.0, 0.0, lc};
//Point on the North pole
Point(2) = {0.0, 0.0, radius, lc};
//Point on the South pole
Point(3) = {0.0, 0.0, -radius, lc};
//Point at 0 deg E, 0 deg N
Point(4) = {0.0, radius, 0.0, lc};
//Point at 180 deg E, 0 deg N
Point(5) = {0.0, -radius, 0.0, lc};
//Point at 90 deg E, 0 deg N
Point(6) = {radius, 0.0, 0.0, lc};
//Point at 90 deg W, 90 deg N
Point(7) = {-radius, 0.0, 0.0, lc};

Circle(1) = {2, 1, 5};
Circle(2) = {5, 1, 6};
Circle(3) = {6, 1, 2};
Circle(4) = {2, 1, 7};
Circle(5) = {5, 1, 7};
Circle(6) = {6, 1, 3};
Circle(7) = {6, 1, 4};
Circle(8) = {4, 1, 3};
Circle(9) = {2, 1, 4};
Circle(10) = {7, 1, 4};
Circle(11) = {7, 1, 3};
Circle(12) = {3, 1, 5};

Line Loop(13) = {1, 2, 3};
Ruled Surface(14) = {13};
Line Loop(15) = {3, 9, -7};
Ruled Surface(16) = {15};
Line Loop(17) = {9, -10, -4};
Ruled Surface(18) = {17};
Line Loop(19) = {4, -5, -1};
Ruled Surface(20) = {19};
Line Loop(21) = {12, 2, 6};
Ruled Surface(22) = {21};
Line Loop(23) = {6, -8, -7};
Ruled Surface(24) = {23};
Line Loop(25) = {8, -11, 10};
Ruled Surface(26) = {25};
Line Loop(27) = {11, 12, 5};
Ruled Surface(28) = {27};
Physical Surface(29) = {16, 14, 20, 18, 28, 26, 24, 22};

//Point on the ocean floor of the North pole
Point(8) = {0.0, 0.0, radius-400000, lc};
//Point on the ocean floor of the South pole
Point(9) = {0.0, 0.0, -(radius-400000), lc};
//Point on the ocean floor at 0 deg E, 0 deg N
Point(10) = {0.0, radius-400000, 0.0, lc};
//Point on the ocean floor at 180 deg E, 0 deg N
Point(11) = {0.0, -(radius-400000), 0.0, lc};
//Point on the ocean floor at 90 deg E, 0 deg N
Point(12) = {radius-400000, 0.0, 0.0, lc};
//Point on the ocean floor at 90 deg W, 90 deg N
Point(13) = {-(radius-400000), 0.0, 0.0, lc};

Circle(31) = {11, 1, 12};
Circle(32) = {12, 1, 10};
Circle(33) = {10, 1, 13};
Circle(34) = {13, 1, 11};
Circle(35) = {11, 1, 9};
Circle(36) = {9, 1, 10};
Circle(37) = {10, 1, 8};
Circle(38) = {8, 1, 11};
Circle(39) = {12, 1, 8};
Circle(40) = {8, 1, 13};
Circle(41) = {13, 1, 9};
Circle(42) = {9, 1, 12};
Line Loop(43) = {32, 37, -39};
Ruled Surface(44) = {43};
Line Loop(45) = {39, 38, 31};
Ruled Surface(46) = {45};
Line Loop(47) = {40, 34, -38};
Ruled Surface(48) = {47};
Line Loop(49) = {40, -33, 37};
Ruled Surface(50) = {49};
Line Loop(51) = {32, -36, 42};
Ruled Surface(52) = {51};
Line Loop(53) = {33, 41, 36};
Ruled Surface(54) = {53};
Line Loop(55) = {34, 35, -41};
Ruled Surface(56) = {55};
Line Loop(57) = {42, -31, 35};
Ruled Surface(58) = {57};
Physical Surface(59) = {48, 46, 44, 50, 54, 56, 58, 52};


Surface Loop(60) = {16, 14, 20, 18, 26, 24, 22, 28};
Surface Loop(61) = {44, 52, 54, 50, 48, 56, 58, 46};
Volume(62) = {60, 61};
Physical Volume(63) = {62};
