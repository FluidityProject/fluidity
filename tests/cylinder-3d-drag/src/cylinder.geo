Point (1) = {0.2, 0.2, 0, 0.005};
Point (2) = {0.25, 0.2, 0, 0.005};
Point (3) = {0.2, 0.25, 0, 0.005};
Point (4) = {0.15, 0.2, 0, 0.005};
Point (5) = {0.2, 0.15, 0, 0.005};
Point (6) = {0, 0, 0, 0.01};
Point (7) = {2.2, 0, 0, 0.05};
Point (8) = {2.2, 0.41, 0, 0.05};
Point (9) = {0, 0.41, 0, 0.01};
Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line (5) = {9, 6};
Line (6) = {9, 8};
Line (7) = {6, 7};
Line (8) = {8, 7};
Line Loop (10) = {6, 8, -7, -5};
Line Loop(101) = {1,2,3,4};
Ruled Surface(102) = {10,101};
Extrude {0,0,0.41} {
  Surface{102};
}
Physical Surface(145) = {144,123,102,115};
Physical Surface(146) = {119};
Physical Surface(147) = {127};
Physical Surface(148) = {131,139,143,135};
Physical Volume(149) = {1};
