SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.5;
Mesh.CharacteristicLengthMax = 0.5;

L = 10;
Point(1) = {-L/2, -L/2, 0};
Point(2) = {L/2, -L/2, 0};
Point(3) = {L/2, L/2, 0};
Point(4) = {-L/2, L/2, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Periodic Line{1} = {-3};
Periodic Line{2} = {-4};

r = 3.0;
Disk(2) = {0, 0, 0, r};
Disk(3) = {-L/2, -L/2, 0, r};
Disk(4) = {L/2, -L/2, 0, r};
Disk(5) = {L/2, L/2, 0, r};
Disk(6) = {-L/2, L/2, 0, r};

f1() = BooleanDifference{ Surface{1}; }{ Surface{2:6}; };
f2() = BooleanIntersection{ Surface{1}; Delete; }{ Surface{2:6}; Delete; };
BooleanUnion{ Surface{ f1() }; Delete; }{ Surface{ f2() }; Delete; }

//Physical Surface("0") = {1};
//Physical Surface("1") = {4, 3, 6, 5, 2};
