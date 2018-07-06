ld=0.5;
lr=2*N; // number of cells in the radial direction
lti=3*N; // number of cells in a pi/4 segment in the tangential direction of the inner radius
lto=lti; // on the outer radius, has to be the same as on the inside for a transfinite surface
ltl=lti; // on the layer, has to be the same for a transfinite surface
p=1;
// Bottom layer
Point(1) = {1.22, 0, 0, 1.0};
Extrude {1.-ld, 0, 0} {
  Point{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{6};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{10};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{14};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{18};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{22};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{26};
}
// Top layer
Extrude {0.5, 0, 0} {
  Point{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{34};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{35};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{39};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{43};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{47};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{51};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{55};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{59};
}
Transfinite Line{1} = lr+1 Using Bump p;
Transfinite Line{2} = lr+1 Using Bump p;
Transfinite Line{6} = lr+1 Using Bump p;
Transfinite Line{10} = lr+1 Using Bump p;
Transfinite Line{14} = lr+1 Using Bump p;
Transfinite Line{18} = lr+1 Using Bump p;
Transfinite Line{22} = lr+1 Using Bump p;
Transfinite Line{26} = lr+1 Using Bump p;
Transfinite Line{3} = lti+1;
Transfinite Line{7} = lti+1;
Transfinite Line{11} = lti+1;
Transfinite Line{15} = lti+1;
Transfinite Line{19} = lti+1;
Transfinite Line{23} = lti+1;
Transfinite Line{27} = lti+1;
Transfinite Line{31} = lti+1;
Transfinite Line{4} = ltl+1;
Transfinite Line{8} = ltl+1;
Transfinite Line{12} = ltl+1;
Transfinite Line{16} = ltl+1;
Transfinite Line{20} = ltl+1;
Transfinite Line{24} = ltl+1;
Transfinite Line{28} = ltl+1;
Transfinite Line{32} = ltl+1;
Transfinite Surface{5} = {1,2,3,4} Alternate;
Transfinite Surface{9} = {3,4,9,10} Alternate;
Transfinite Surface{13} = {9,10,11,12} Alternate;
Transfinite Surface{17} = {11,12,13,14} Alternate;
Transfinite Surface{21} = {13,14,15,16} Alternate;
Transfinite Surface{25} = {15,16,17,18} Alternate;
Transfinite Surface{29} = {17,18,19,20} Alternate;
Transfinite Surface{33} = {19,20,1,2} Alternate;
Transfinite Surface{38} = {2, 21, 4, 23} Alternate;
Transfinite Surface{42} = {4, 23, 10, 25} Alternate;
Transfinite Surface{46} = {10,25, 12, 27} Alternate;
Transfinite Surface{50} = {12, 27, 14, 29} Alternate;
Transfinite Surface{54} = {14, 29, 16, 31} Alternate;
Transfinite Surface{58} = {16, 31, 18, 33} Alternate;
Transfinite Surface{62} = {18, 33, 20, 35} Alternate;
Transfinite Surface{66} = {20, 35, 2, 21} Alternate;
Transfinite Line{37, 41, 45, 49, 53, 57, 61, 65} = lto+1;
Transfinite Line{34, 35, 39, 43, 47, 51, 55, 59} = lr+1 Using Bump p;

// inner radius
Physical Line(1) = {3, 7, 11, 15, 19, 23, 27, 31};
// outer radius
Physical Line(2) = {37, 41, 45, 49, 53, 57, 61, 65};
// layer
Physical Line(3) = {4, 8, 12, 16, 20, 24, 28, 32};

// "upper mantle"
Physical Surface(1) = {38, 42, 46, 50, 54, 58, 62, 66};
// "lower mantle"
Physical Surface(2) = {5, 9, 13, 17, 21, 25, 29, 33};
