lc = 1.5;
L = 10; // length of the square
R = 5.0; // circle radius
TFpoints = 5;
//+
Point(1) = {-L, -L, 0, lc};
//+
Point(2) = {-L, L, 0, lc};
//+
Point(3) = {L, -L, 0, lc};
//+
Point(4) = {L, L, 0, lc};
//+
Point(5) = {0, 0, 0, lc};
//+
Point(6) = {R, 0, 0, lc};
//+
Point(7) = {0., R, 0, lc};
//+
Point(8) = {0., -R, 0, lc};
//+
Point(9) = {-R, 0., 0, lc};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 2};
//+
Circle(5) = {9, 5, 7};
//+
Circle(6) = {7, 5, 6};
//+
Circle(7) = {6, 5, 8};
//+
Circle(8) = {8, 5, 9};
//+
Line Loop(9) = {5, 6, 7, 8};
//+
Plane Surface(0) = {9};
//+
Line Loop(11) = {1, 2, 3, 4};
//+
Plane Surface(1) = {9, 11};

Physical Surface(0) = {0};
//+
Physical Surface(1) = {1};

//+
Transfinite Line {1, 3} = TFpoints Using Progression 1;
