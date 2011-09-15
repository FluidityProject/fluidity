// In xy-plane, z unused
// Define km
km=1000.0;

aspect=1;
//aspect=1.0;
// Ice shelf front depth
iceshelffrontdepth=0.2*km;
// Ice shelf depth below front (draft + additional depth below max draft)
oceandepthminusfront=0.9*km;
// Shelf length
shelflength=0.55*km * aspect;
// Open ocean length
oceanlength=0.55*km * aspect;
// Inlet and outlet height
inoutletlength=0.3*km;

oceandeptha=0.3*km;
oceandepthb=0.5*km;
//oceandepthc=0.1*km + icesheffrontdepth;
oceandepthc=0.3*km;


// Southern boundary
Point(1) = {0,0,0};
Extrude {0, 0, oceandeptha} { Point{1}; Layers{6}; }
Extrude {0, 0, oceandepthb} { Point{2}; Layers{10}; }
Extrude {0, 0, oceandepthc} { Point{3}; Layers{6}; }

width=1.1*km;
//ylayers=22;
ylayers=4;



// // Extrude north in latitude
// //Line(1) = {1,2};
Extrude {shelflength,0,0} { Line{1}; Line{2}; Line{3}; Layers{20}; }
// 
// // Extrude north in latitude again
// //Line(2) = {3,4};
// //Line(3) = {1,2};
// //Line(4) = {2,3};
Extrude {oceanlength,0,0} { Line{4}; Line{8}; Line{12}; Layers{20}; }
// 
// // Extrude up vertically, creating ice shelf front
// //Line(6) = {5,6};
// //Line(7) = {3,5};
// //Line(8) = {4,6};
//Extrude {0, 0, iceshelffrontdepth} { Line{26}; Layers{4}; }
// 
// //Line(10) = {4,7};
// //Line(11) = {6,8};
// //Line(12) = {7,8};
// 
// 


Extrude {0, width, 0} {
  Surface{11, 23, 7, 19, 15, 27, 31}; Layers{ylayers};
}

// // Boundaries
// // Ocean free surface
// Physical Line(101) = {28};
// // Bottom
// Physical Line(102) = {5,17};
// // South
// Physical Line(103) = {1,2,3};
// 
// // North
// Physical Line(104) = {20};
// // North inlet
// Physical Line(1041) = {16};
// // North outlet
// Physical Line(1042) = {24,30};
// 
// // Ice front
// Physical Line(105) = {29};
// // Ice slope
// Physical Line(106) = {14};
// 
// Physical Surface(107) = {7,11,15,19,23,27,31};
// 



// Boundaries
// Ocean free surface
Physical Surface(101) = {150};
// Bottom
Physical Surface(102) = {92,114};
// South
Physical Surface(103) = {80,36,124};

// North
//Physical Surface(104) = {70};
// North inlet
//Physical Surface(1041) = {114};
// North outlet
//Physical Surface(1042) = {158,176};

Physical Surface(104) = {66,110,154};

// Ice front
//Physical Surface(105) = {184};
// Ice slope
Physical Surface(106) = {128};

// West
Physical Surface(107) = {49,71,93,115,137,159};
// East
Physical Surface(108) = {7,11,15,19,23,27};

// Volume
Physical Volume(110) = {1,2,3,4,5,6};

