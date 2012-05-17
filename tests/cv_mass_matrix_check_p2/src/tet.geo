/* Simple tetrahedral domain with one element
*/

// define variable used for characteristic length
lc = 2.0;

// define geometric entities

// define points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {0.2, 0.8, 0.0, lc};
Point(4) = {0.3, 0.8, 2.0, lc};

//define lines
Line (1) = {1, 2};
Line (2) = {1, 3};
Line (3) = {1, 4};
Line (4) = {2, 3};
Line (5) = {2, 4};
Line (6) = {3, 4};

//define line loop
Line Loop (8)  = {2,  6, -3};
Line Loop (10) = {2, -4, -1};
Line Loop (12) = {3, -5, -1};
Line Loop (14) = {4,  6, -5};

// define surface
Plane Surface (8)  = {8};
Plane Surface (10) = {10};
Plane Surface (12) = {12};
Plane Surface (14) = {14};

// define surface loop
Surface Loop (16) = {14, 10, 8, 12};

// define volume
Volume (16) = {16};

// define physical entities

// define physicial surface
Physical Surface (101) = {8, 10, 12, 14};

// define physical volume
Physical Volume (102) = {16};
