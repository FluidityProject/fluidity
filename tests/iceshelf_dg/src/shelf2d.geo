// In xy-plane, z unused
// Define km
km=1000.0;

//aspect=1000.0;
aspect=1.0;
// Ice shelf front depth
iceshelffrontdepth=0.2*km;
// Ice shelf depth below front (draft + additional depth below max draft)
oceandepthminusfront=0.9*km;
// Shelf length
shelflength=0.55*km * aspect;
// Open ocean length
oceanlength=0.55*km * aspect;


// Southern boundary
Point(1) = {0,0,0};
Extrude {0, oceandepthminusfront,0} {
Point{1}; Layers{36};
}
// Extrude north in latitude
//Line(1) = {1,2};
Extrude {shelflength,0,0} {
Line{1}; Layers{50};
}
// Extrude north in latitude again
//Line(2) = {3,4};
//Line(3) = {1,2};
//Line(4) = {2,3};
Extrude {oceanlength,0,0} {
Line{2}; Layers{50};
}
// Extrude up vertically, creating ice shelf front
//Line(6) = {5,6};
//Line(7) = {3,5};
//Line(8) = {4,6};
Extrude {0, iceshelffrontdepth,0} {
Line{8}; Layers{8};
}
//Line(10) = {4,7};
//Line(11) = {6,8};
//Line(12) = {7,8};

// Boundaries
// Ocean free surface
Physical Line(101) = {10};
// Bottom
Physical Line(102) = {3,7};
// South
Physical Line(103) = {1};
// North
Physical Line(104) = {6,12};
// Ice front
Physical Line(105) = {11};
// Ice slope
Physical Line(106) = {4};

Physical Surface(107) = {13,9,5};

