lc = 2;
R = 7;
L = 10;
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {R,0.0,0.0,lc};
Point(3) = {0,R,0.0,lc};
Circle(1) = {2,1,3};
Point(4) = {-R,0,0.0,lc};
Point(5) = {0,-R,0.0,lc};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Point(6) = {0,0,-R,lc};
Point(7) = {0,0,R,lc};
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};
Line Loop(13) = {-2,-8,10};
Ruled Surface(14) = {13};
Line Loop(15) = {10,3,7};
Ruled Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Ruled Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Ruled Surface(20) = {19};
Line Loop(21) = {5,12,1};
Ruled Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Ruled Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Ruled Surface(26) = {25};
Line Loop(27) = {4,-12,6};
Ruled Surface(28) = {27};

Surface Loop(29) = {-28,26,16,-14,20,24,-22,18};

Point(8) = {L, L, -L, lc};
//+
Point(9) = {-L, L, -L, lc};
//+
Point(10) = {-L, -L, -L, lc};
//+
Point(11) = {L, -L, -L, lc};
//+
Line(31) = {8, 11};
//+
Line(32) = {11, 10};
//+
Line(33) = {10, 9};
//+
Line(34) = {9, 8};

//+
Point(12) = {L, L, L, lc};
//+
Point(13) = {-L, L, L, lc};
//+
Point(14) = {L, -L, L, lc};
//+
Point(15) = {-L, -L, L, lc};
//+
Line(35) = {13, 12};
//+
Line(36) = {12, 14};
//+
Line(37) = {14, 15};
//+
Line(38) = {15, 13};
//+
Line(39) = {9, 13};
//+
Line(40) = {10, 15};
//+
Line(41) = {8, 12};
//+
Line(42) = {11, 14};
//+
Line Loop(43) = {34, 41, -35, -39};
//+
Plane Surface(44) = {43};
//+
Line Loop(45) = {33, 39, -38, -40};
//+
Plane Surface(46) = {45};
//+
Line Loop(47) = {32, 33, 34, 31};
//+
Plane Surface(48) = {47};
//+
Line Loop(49) = {31, 42, -36, -41};
//+
Plane Surface(50) = {49};
//+
Line Loop(51) = {32, 40, -37, -42};
//+
Plane Surface(52) = {51};
//+
Line Loop(53) = {37, 38, 35, 36};
//+
Plane Surface(54) = {53};
//+
Surface Loop(55) = {48, 52, 46, 44, 50, 54};
//+
Volume(0) = {29, 55};
//+
Volume(1) = {29};
//+
Transfinite Line {32, 34, 35, 37} = 10 Using Progression 1;
//+
Transfinite Line {33, 31, 36, 38} = 10 Using Progression 1;
//+
Transfinite Surface {48};
//+
Transfinite Surface {54};
//+
Transfinite Surface {50};
//+
Transfinite Surface {46};
//+
Transfinite Surface {52};
//+
Transfinite Surface {44};
