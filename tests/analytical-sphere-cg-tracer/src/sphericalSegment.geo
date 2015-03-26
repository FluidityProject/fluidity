Printf("Reading function and point definitions...");

Function DrawParallel
 alpha = parallelSectionStartingY/parallelSectionStartingX;
 parallelRadius = Sqrt(parallelSectionStartingX^2 + parallelSectionStartingY^2);
 deltaAlpha = (parallelSectionEndingY/parallelSectionEndingX - parallelSectionStartingY/parallelSectionStartingX)/pointsOnParallel;
 For t In {1:pointsOnParallel}
  alpha += deltaAlpha;
  newParallelPointX = (parallelSectionStartingX/Fabs(parallelSectionStartingX))*parallelRadius/(Sqrt(alpha^2 + 1));
  newParallelPointY = newParallelPointX*alpha;
  newPointOnParallel = newp; Point(newPointOnParallel) = {newParallelPointX, newParallelPointY, 0};
 EndFor
 BSpline(newParallelID) = {firstPointOnParallel, newPointOnParallel-(pointsOnParallel-1):newPointOnParallel, lastPointOnParallel};
Return

//The coordinates of various points on the Surface of the Earth are converted to radians
// and stereographic coordinates below. The geometry drawn is a square (in lon-lat space)
// centered about the intersection of the Equator and the Greenwich meridian.

Point_longitude_rad = (01 + (30/60))*(Pi/180);
Point_latitude_rad =  (01 + (30/60))*(Pi/180);
cell_layers = 5; //Number of element-layers per direction.

earth_surface_radius = 6.37101e+06;

Printf("Defining points...");
//Point A is located at the North-East corner of the domain
Point_A_longitude_rad = Point_longitude_rad;
Point_A_latitude_rad = Point_latitude_rad;
Point_A_stereographic_x = -Cos(Point_A_longitude_rad)*Cos(Point_A_latitude_rad)/( 1 + Sin(Point_A_latitude_rad) );
Point_A_stereographic_y = -Cos(Point_A_latitude_rad)*Sin(Point_A_longitude_rad)/( 1 + Sin(Point_A_latitude_rad) );

//Point B is located at the North-West corner of the domain
Point_B_longitude_rad = -Point_longitude_rad;
Point_B_latitude_rad = Point_latitude_rad;
Point_B_stereographic_x = -Cos(Point_B_longitude_rad)*Cos(Point_B_latitude_rad)/( 1 + Sin(Point_B_latitude_rad) );
Point_B_stereographic_y = -Cos(Point_B_latitude_rad)*Sin(Point_B_longitude_rad)/( 1 + Sin(Point_B_latitude_rad) );

//Point C is located at the South-West corner of the domain
Point_C_longitude_rad = -Point_longitude_rad;
Point_C_latitude_rad = -Point_latitude_rad;
Point_C_stereographic_x = -Cos(Point_C_longitude_rad)*Cos(Point_C_latitude_rad)/( 1 + Sin(Point_C_latitude_rad) );
Point_C_stereographic_y = -Cos(Point_C_latitude_rad)*Sin(Point_C_longitude_rad)/( 1 + Sin(Point_C_latitude_rad) );

//Point D is located at the South-East corner of the domain
Point_D_longitude_rad = Point_longitude_rad;
Point_D_latitude_rad = -Point_latitude_rad;
Point_D_stereographic_x = -Cos(Point_D_longitude_rad)*Cos(Point_D_latitude_rad)/( 1 + Sin(Point_D_latitude_rad) );
Point_D_stereographic_y = -Cos(Point_D_latitude_rad)*Sin(Point_D_longitude_rad)/( 1 + Sin(Point_D_latitude_rad) );

//Distance calculation, along loxodrome of constant latitude (lines AB & CD):
// see also http://www.movable-type.co.uk/scripts/latlong.html
delta = Log(Tan(Point_B_latitude_rad/2 + Pi/4)/(Tan(Point_A_latitude_rad/2 + Pi/4)));
q = Cos(Point_A_latitude_rad);
distance = q*(Point_A_longitude_rad - Point_B_longitude_rad)*earth_surface_radius;
elementSize = distance/cell_layers;
Printf("Characteristic size of the domain: %g m", distance);
Printf("Expected element edge length: %g m", elementSize);

//Create point, line, line-loop surface and field ID's
IP = newp;
IL = newl;
ILL = newll;
IS = news;
IFI = newf;

//Draw a point on the centre of the earth and one on the North pole, then create a sphere using
// these points, note that remaining points will be drawin in the stereographic projection plane.
Point ( IP + 0 ) = {0, 0, 0 };
Point ( IP + 1 ) = {0, 0, earth_surface_radius};
PolarSphere ( IS + 0 ) = {IP , IP+1};

Printf("Creating the domain...");
Point(IP + 60000) = {Point_A_stereographic_x, Point_A_stereographic_y, 0};
Point(IP + 60001) = {Point_B_stereographic_x, Point_B_stereographic_y, 0};
Point(IP + 60002) = {Point_C_stereographic_x, Point_C_stereographic_y, 0};
Point(IP + 60003) = {Point_D_stereographic_x, Point_D_stereographic_y, 0};

BSpline(IL + 1001) = {IP + 60000, IP + 60003};

//Draw North-most parallel of the domain, through points A, B
pointsOnParallel = 200; //Assign parameters to variables, then call function DrawParallel,
parallelSectionStartingX = Point_A_stereographic_x;
parallelSectionStartingY = Point_A_stereographic_y;
firstPointOnParallel = IP + 60000;
parallelSectionEndingX = Point_B_stereographic_x;
parallelSectionEndingY = Point_B_stereographic_y;
lastPointOnParallel = IP + 60001;
newParallelID = IL + 1002;
Call DrawParallel;

//Draw North-most parallel of the domain, through points C, D
pointsOnParallel = 200; //Assign parameters to variables, then call function DrawParallel,
parallelSectionStartingX = Point_D_stereographic_x;
parallelSectionStartingY = Point_D_stereographic_y;
firstPointOnParallel = IP + 60003;
parallelSectionEndingX = Point_C_stereographic_x;
parallelSectionEndingY = Point_C_stereographic_y;
lastPointOnParallel = IP + 60002;
newParallelID = IL + 1003;
Call DrawParallel;

BSpline(IL + 1004) = {IP + 60001, IP + 60002};
Transfinite Line{IL + 1001, IL + 1002, IL + 1003, IL + 1004} = cell_layers + 1;

Line Loop(1005) = {1002, 1004, -1005, -1003};
Ruled Surface (1006) = {1005} In Sphere {IS + 0};
Transfinite Surface{1006};

Printf("Declaring surface on sphere to be meshed...");
Physical Surface(1008) = {1006};

Physical Line(1010) = {1002}; //Physical line corresponding to meridian segment through points A & D.
Physical Line(1011) = {1005}; //Physical line corresponding to meridian segment through points B & C.
Physical Line(1012) = {1003}; //Physical line corresponding to North-most parallel.
Physical Line(1009) = {1004}; //Physical line corresponding to Sorth-most parallel.

//Printf("Assigning characteristic mesh sizes...");
//Field[ IFI + 1] = MathEval;
//Field[ IFI + 1].F = Sprintf("%g", elementSize);

//Background Field = IFI + 1;

// Dont extent the elements sizes from the boundary inside the domain
//Mesh.CharacteristicLengthExtendFromBoundary = 0;

//Set some options for better png output
General.Color.Background = {255,255,255};
General.Color.BackgroundGradient = {255,255,255};
General.Color.Foreground = Black;
Mesh.Color.Lines = {0,0,0};
