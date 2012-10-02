//

lf = 0.001;
lm = 0.015;
lc = 0.075;
mc = 0.0075;


// box points:
Point (6) = {0, 0, 0, lc};
Point (7) = {2.2, 0, 0, lc};
Point (8) = {2.2, 0.41, 0, lc};
Point (9) = {0, 0.41, 0, lc};
Line (5) = {9, 6};
Line (6) = {9, 8};
Line (7) = {6, 7};
Line (8) = {8, 7};
// Box surface:
Line Loop(9) = {6, 8, -7, -5};
Plane Surface(10) = {9};



// Physical IDs:
//front/inlet:
Physical Line(1) = {5};
//top and bottom
Physical Line(2) = {6, 7};
//back/outlet:
Physical Line(3) = {8};
// Surface
Physical Surface(11) = {10};


