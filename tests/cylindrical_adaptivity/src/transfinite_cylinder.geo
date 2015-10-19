lr=2; // number of cells in the radial direction
lti=lr; // number of cells in a pi/4 segment in the tangential direction of the inner radius
lto=lti; // has to be for a transfinite surface
p=1;
Point(1) = {1.22, 0, 0, 1.0};
Extrude {1, 0, 0} {
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
Transfinite Line{4} = lto+1;
Transfinite Line{8} = lto+1;
Transfinite Line{12} = lto+1;
Transfinite Line{16} = lto+1;
Transfinite Line{20} = lto+1;
Transfinite Line{24} = lto+1;
Transfinite Line{28} = lto+1;
Transfinite Line{32} = lto+1;
Transfinite Surface{5} = {1,2,3,4} Alternate;
Transfinite Surface{9} = {3,4,9,10} Alternate;
Transfinite Surface{13} = {9,10,11,12} Alternate;
Transfinite Surface{17} = {11,12,13,14} Alternate;
Transfinite Surface{21} = {13,14,15,16} Alternate;
Transfinite Surface{25} = {15,16,17,18} Alternate;
Transfinite Surface{29} = {17,18,19,20} Alternate;
Transfinite Surface{33} = {19,20,1,2} Alternate;
// inner radius
Physical Line(1) = {3, 7, 11, 15, 19, 23, 27, 31};
// outer radius
Physical Line(2) = {4, 8, 12, 16, 20, 24, 28, 32};
Physical Surface(0) = {5, 9, 13, 17, 21, 25, 29, 33};
