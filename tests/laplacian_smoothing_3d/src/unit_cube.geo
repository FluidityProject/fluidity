char_len = 0.1;
//This produces F4 and serves as the \Omega_{c} in the 2D experiments
Point(1) = {0, 0, 0, char_len};
Point(2) = {1, 0, 0, char_len};
Point(3) = {1, 1, 0, char_len};
Point(4) = {0, 1, 0, char_len};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

//Produce F2
Point(5) = {0, 0, 1, char_len};
Point(6) = {1, 0, 1, char_len};
Point(7) = {1, 1, 1, char_len};
Point(8) = {0, 1, 1, char_len};
Line(5) = {8, 7};
Line(6) = {7, 6};
Line(7) = {6, 5};
Line(8) = {5, 8};
Line Loop(7) = {5, 6, 7, 8};
Plane Surface(8) = {7};
//Produce F1
Line(9) = {4,8};
Line(10) = {5,1};
Line Loop(9) = {9, -8, 10, 4};
Plane Surface(10) = {9};
//Produce F3
Line(11) = {3, 7};
Line(12) = {6, 2};
Line Loop(13) = {11, 6, 12, -2};
Plane Surface(14) = {13};
//Produce F5
Line Loop(15) = {-12, 7, 10, -3};
Plane Surface(16) = {15};
//Produce F6
Line Loop(17) = {1, 11, -5, -9};
Plane Surface (18) = {17};

//Master Surface Loop
Surface Loop(19) = {6, 8, 10, 14, 16, 18};
Volume(20) = {19};


//Labelling of Points to make processing easier for C and Python code
Physical Point(21) = {1,2,3,4,5,6,7,8};
Physical Line(22) = {1,2,3,4,5,6,7,8};
Physical Surface(23) = {6,8,10,14,16,18};
Physical Volume(24) = {20};
