IP = newp;
IL = newl;
IS = news;IF = newf;
Point(IP + 0) = {0, 0, 0};
Point(IP + 1) = {0, 0, 6.371010e+06};
PolarSphere(IS) = {IP + 0, IP + 1};

// Four point, making a rectangle 10x5 degrees in size over north pole
// NP: 0,90
// P1: -2.5,85
// P2: 2.5,85
// P3: -180+2.5,85
// P4: 180-2.5,85

Point(IP + 2) = {-0.043619387365335986, 0.0019044635814622068, 0.};
Point(IP + 3) = {-0.043619387365335986, -0.0019044635814622016, 0.};
Point(IP + 4) = {0.043619387365335986, 0.0019044635814622068, 0.};
Point(IP + 5) = {0.043619387365335986, -0.0019044635814622125, 0.};

Line(1) = {3, 4};
Line(2) = {4, 6};
Line(3) = {6, 5};
Line(4) = {5, 3};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Line(7) = {1, 2, 3, 4};
Physical Surface(8) = {6};
Field[1] = MathEval;
Background Field = 1;
Field[1].F = "10000";
