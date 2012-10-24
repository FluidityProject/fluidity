/* Simple triangle domain with one element
*/

// define variable used for characteristic length
lc = 1.0;

// define geometric entities

// define points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {0.2, 0.8, 0.0, lc};

//define lines
Line(4) = {1,2};
Line(5) = {2,3};
Line(6) = {3,1};

//define line loop
Line Loop(7) = {4,5,6};

// define surface
Plane Surface(8) = {7};

// define physical entities

// define physical lines
Physical Line (101) = {4,5,6};

// define physicial surface
Physical Surface(102) = {8};
