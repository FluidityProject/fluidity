Point (1) = {0, 0.25, 0, 0.05};
Point (2) = {0.5, 0.25, 0, 0.05};
Point (3) = {0.5, 0.75, 0, 0.05};
Point (4) = {0, 0.75, 0, 0.05};
Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 1};
Line Loop (5) = {1, 2, 3, 4};
Plane Surface (1) = {5};
Physical Line (1) = {1};
Physical Line (2) = {2};
Physical Line (3) = {3};
Physical Line (4) = {4};
Extrude {0,0,0.05} {
  Surface{1}; Layers{1};
}
Physical Surface(29) = {14};
Physical Surface(30) = {18};
Physical Surface(31) = {22};
Physical Surface(32) = {26};
Physical Surface(33) = {27};
Physical Surface(34) = {1};
Physical Volume(1) = {1};
