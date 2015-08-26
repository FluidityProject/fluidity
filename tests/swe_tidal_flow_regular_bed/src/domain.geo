Lx = 14000.0;
nx = 1000; // Number of points in the x direction.
dx = Lx/nx;
Point(1) = {0,0,0,dx};
Extrude {Lx,0,0} {
  Point{1}; Layers{nx};
}
Extrude {0,1000,0} {
  Line{1}; Layers{1};
}

// Reserve 1 and 2 for top and bottom of extruded mesh.

// Outer ends of the channel (x=0 and x=Lx)
Physical Line(3) = {3};
Physical Line(4) = {4};

// Sides of the channel (y=0 and y=1e3)
Physical Line(5) = {1};
Physical Line(6) = {2};
 
Physical Surface(1) = {5};

