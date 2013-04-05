earth_surface_radius = 6.37101e+06;
elementSize_deg = 10;
elementSize_rad = elementSize_deg*(Pi/180);
elementSize_m = 2*elementSize_rad*earth_surface_radius;
Printf("Expected element edge length: %g deg", elementSize_deg);
Printf("Expected element edge length: %g rad", elementSize_rad);
Printf("Expected element edge length: %g m", elementSize_m);

//Create point, line, line-loop surface and field ID's
IP = newp;
IFI = newf;

//Draw a point on the centre of the earth and one on the North pole, then create a sphere using
// these points, note that remaining points will be drawin in the stereographic projection plane.
Point ( IP + 0 ) = {0, 0, 0 };
Point ( IP + 1 ) = {0, 0, earth_surface_radius};
Point ( IP + 2 ) = {earth_surface_radius, 0, 0};
Point ( IP + 3 ) = {0, 0, -earth_surface_radius};


Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{1, 2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{3, 6};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{9, 12};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{15, 18};
}

Physical Surface(1008) = {8, 5, 11, 14, 17, 20, 23, 26};

//Printf("Assigning characteristic mesh sizes...");
Field[ IFI + 1] = MathEval;
Field[ IFI + 1].F = Sprintf("%g", elementSize_m);

Background Field = IFI + 1;

// Dont extent the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;

//Set some options for better png output
General.Color.Background = {255,255,255};
General.Color.BackgroundGradient = {255,255,255};
General.Color.Foreground = Black;
Mesh.Color.Lines = {0,0,0};
