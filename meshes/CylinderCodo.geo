// Gmsh project created on Mon Jan 29 14:40:11 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, -0.1, 1.0};
//+
Point(2) = {1, 0, -0.1, 1.0};
//+
Point(3) = {1, 1, -0.1, 1.0};
//+
Point(4) = {0, 1, -0.1, 1.0};
//+
Point(5) = {0, 0, 1.1, 1.0};
//+
Point(6) = {1, 0, 1.1, 1.0};
//+
Point(7) = {1, 1, 1.1, 1.0};
//+
Point(8) = {0, 1, 1.1, 1.0};
//+
Circle(15) = {0, 0.5, 0.75, 0.1, 0, 2*Pi};
//+
Circle(16) = {1, 0.5, 0.25, 0.1, 0, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0.5, 0.75}, Pi/2} {
  Curve{15}; 
}
//+
Rotate {{0, 1, 0}, {1, 0.5, 0.25}, Pi/2} {
  Curve{16}; 
}
//+
Line(17) = {4, 1};
//+
Line(18) = {1, 5};
//+
Line(19) = {5, 8};
//+
Line(20) = {8, 4};
//+
Line(21) = {1, 2};
//+
Line(22) = {2, 6};
//+
Line(23) = {6, 7};
//+
Line(24) = {7, 3};
//+
Line(25) = {3, 2};
//+
Line(26) = {4, 3};
//+
Line(27) = {6, 5};
//+
Line(28) = {7, 8};
//+

//+
Curve Loop(1) = {16};
//+
Curve Loop(2) = {22, 23, 24, 25};
//+
Curve Loop(3) = {16};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {15};
//+
Plane Surface(1) = {4};
//+
Curve Loop(5) = {16};
//+
Curve Loop(6) = {25, 22, 23, 24};
//+
Plane Surface(3) = {5, 6};
//+
Curve Loop(7) = {17, 18, 19, 20};
//+
Curve Loop(8) = {15};
//+
Plane Surface(11) = {7, 8};
//+
Curve Loop(9) = {24, -26, -20, -28};
//+
Plane Surface(5) = {9};
//+
Curve Loop(10) = {21, 22, 27, -18};
//+
Plane Surface(6) = {10};
//+
Curve Loop(11) = {25, -21, -17, 26};
//+
Plane Surface(7) = {11};
//+
Curve Loop(12) = {28, -19, -27, 23};
//+
Plane Surface(8) = {12};
//+
Surface Loop(1) = {1, 3, 7, 6, 8, 5, 11, 2};
//+
Volume(1) = {1};
