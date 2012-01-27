// define a characteristic length
lc = 60.0;

len = 300.0;

// define points
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {len,0.0,0.0,lc};
Point(3) = {len,len,0.0,lc};
Point(4) = {0.0,len,0.0,lc};

// define line's
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// define line loop
Line Loop(1) = {1,2,3,4};

// define surface
Plane Surface(1) = {1};

// Rotate the domain about the z axis by 45 degrees (Pi/4)
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Surface{1};
}

// Tag the lines and surface
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};

Physical Surface(11) = {1};

// Make the grid regular
Transfinite Surface"*";
