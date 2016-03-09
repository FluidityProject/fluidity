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
//and stereographic coordinatesbelow.

//Point A is located 2 deg 00'W, 48 deg 45'N.
Point_A_longitude_rad = -(2 + (00/60))*(Pi/180);
Point_A_latitude_rad = (48 + (45/60))*(Pi/180);
Point_A_stereographic_x = -Cos(Fabs(Point_A_longitude_rad))*Cos(Fabs(Point_A_latitude_rad))/( 1 + Sin(Fabs(Point_A_latitude_rad)) );
Point_A_stereographic_y = Cos(Fabs(Point_A_latitude_rad))*Sin(Fabs(Point_A_longitude_rad))/( 1 + Sin(Fabs(Point_A_latitude_rad)) );

//Point B is located at 9 deg 46'W, 48 deg 45'N
Point_B_longitude_rad = -(9 + (46/60))*(Pi/180);
Point_B_latitude_rad = (48 + (45/60))*(Pi/180);
Point_B_stereographic_x = -Cos(Fabs(Point_B_longitude_rad))*Cos(Fabs(Point_B_latitude_rad))/( 1 + Sin(Fabs(Point_B_latitude_rad)) );
Point_B_stereographic_y = Cos(Fabs(Point_B_latitude_rad))*Sin(Fabs(Point_B_longitude_rad))/( 1 + Sin(Fabs(Point_B_latitude_rad)) );

//Point C is located at 9 deg 45'W, 56 deg 45'N
Point_C_longitude_rad = -(9 + (45/60))*(Pi/180);
Point_C_latitude_rad = (56 + (45/60))*(Pi/180);
Point_C_stereographic_x = -Cos(Fabs(Point_C_longitude_rad))*Cos(Fabs(Point_C_latitude_rad))/( 1 + Sin(Fabs(Point_C_latitude_rad)) );
Point_C_stereographic_y = Cos(Fabs(Point_C_latitude_rad))*Sin(Fabs(Point_C_longitude_rad))/( 1 + Sin(Fabs(Point_C_latitude_rad)) );

//Point D is located at 2 deg 00'W, 56 deg 45'N
Point_D_longitude_rad = -(2 + (00/60))*(Pi/180);
Point_D_latitude_rad = (56 + (45/60))*(Pi/180);
Point_D_stereographic_x = -Cos(Fabs(Point_D_longitude_rad))*Cos(Fabs(Point_D_latitude_rad))/( 1 + Sin(Fabs(Point_D_latitude_rad)) );
Point_D_stereographic_y = Cos(Fabs(Point_D_latitude_rad))*Sin(Fabs(Point_D_longitude_rad))/( 1 + Sin(Fabs(Point_D_latitude_rad)) );

//Create point, line, line-loop surface and field ID's
IP = newp;
IL = newl;
ILL = newll;
IS = news;
IFI = newf;

//Draw a point on the centre of the earth and one on the North pole, then create a sphere using
// these points, note that remaining points will be drawin in the stereographic projection plane.
Point ( IP + 0 ) = {0, 0, 0 };
Point ( IP + 1 ) = {0, 0,6.37101e+06};
PolarSphere ( IS + 0 ) = {IP , IP+1};

Printf("Creating Meridian & Parallel segments to close the domain...");
//So far, only points on the shorelies and curves approximating shorelines have been
//constructed. we now add points in the sea, to close the domain boundaries
Point(IP + 60000) = {Point_A_stereographic_x, Point_A_stereographic_y, 0};
Point(IP + 60001) = {Point_B_stereographic_x, Point_B_stereographic_y, 0};
Point(IP + 60002) = {Point_C_stereographic_x, Point_C_stereographic_y, 0};
Point(IP + 60003) = {Point_D_stereographic_x, Point_D_stereographic_y, 0};

//Draw South-most parallel of the Domain.
pointsOnParallel = 200;               //Assign parameters to variables and then call function DrawParallel,
parallelSectionStartingX = Point_A_stereographic_x;
parallelSectionStartingY = Point_A_stereographic_y;
firstPointOnParallel = IP + 60000;
parallelSectionEndingX = Point_B_stereographic_x;
parallelSectionEndingY = Point_B_stereographic_y;
lastPointOnParallel = IP + 60001;
newParallelID = IL + 1000;
Call DrawParallel;

BSpline(IL + 1001) = {IP + 60001, IP + 60002};
BSpline(IL + 1002) = {IP + 60000, IP + 60003};

////Draw North-most parallel of the Domain.
pointsOnParallel = 200;            //Assign parameters to variables and then call function DrawParallel,
parallelSectionStartingX = Point_C_stereographic_x;
parallelSectionStartingY = Point_C_stereographic_y;
firstPointOnParallel = IP + 60002;
parallelSectionEndingX = Point_D_stereographic_x;
parallelSectionEndingY = Point_D_stereographic_y;
lastPointOnParallel = IP + 60003;
newParallelID = IL + 1003;
Call DrawParallel;

Printf("Declaring surface on sphere to be meshed...");
Line Loop(1006) = {1003, -1004, -1002, -1001};
Plane Surface(1007) = {1006};
Physical Surface(1008) = {1007};

Physical Line(1009) = {1001}; //Physical line corresponding to south-most parallel.
Physical Line(1010) = {1002}; //Physical line corresponding to meridian segment through point A.
Physical Line(1011) = {1003}; //Physical line corresponding to meridian segment through point B.
Physical Line(1013) = {1004}; //Physical line corresponding to north-most parallel.

Printf("Assigning characteristic mesh sizes...");
Field[ IFI + 1] = MathEval;
Field[ IFI + 1].F = "50000";

Background Field = IFI + 1;

// Dont extent the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;

//Set some options for better png output
General.Color.Background = {255,255,255};
General.Color.BackgroundGradient = {255,255,255};
General.Color.Foreground = Black;
Mesh.Color.Lines = {0,0,0};

