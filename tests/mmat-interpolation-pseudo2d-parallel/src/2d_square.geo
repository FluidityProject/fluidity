Point(1) = {0.0,0.0,0,0.062831853071795868};
Point(2) = {3.1415926535897931,0.0,0,0.062831853071795868};
Point(3) = {3.1415926535897931,3.1415926535897931,0,0.062831853071795868};
Point(4) = {0.0,3.1415926535897931,0,0.062831853071795868};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Extrude {0,0,0.02} {
  Surface{6};
}
Physical Surface(8) = {23};
Physical Surface(9) = {19};
Physical Surface(10) = {15};
Physical Surface(11) = {27};
Physical Surface(12) = {28};
Physical Surface(13) = {6};
Physical Volume(35) = {1};
