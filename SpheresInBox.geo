SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.5;
Mesh.CharacteristicLengthMax = 0.5;

r = 4.0;
L = 10;

Box(1) = {0,0,0,L, L, L};

Periodic Surface{2} = {1} Translate{L,0,0};
Periodic Surface{4} = {3} Translate{0,L,0};
Periodic Surface{6} = {5} Translate{0,0,L};

Sphere(2) = {0, 0, 0, r, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {L, L, L, r, -Pi/2, Pi/2, 2*Pi};

Sphere(4) = {0, L, L, r, -Pi/2, Pi/2, 2*Pi};
Sphere(5) = {L, 0, L, r, -Pi/2, Pi/2, 2*Pi};
Sphere(6) = {L, L, 0, r, -Pi/2, Pi/2, 2*Pi};

Sphere(7) = {L, 0, 0, r, -Pi/2, Pi/2, 2*Pi};
Sphere(8) = {0, L, 0, r, -Pi/2, Pi/2, 2*Pi};
Sphere(9) = {0, 0, L, r, -Pi/2, Pi/2, 2*Pi};

Sphere(10) = {L/2, L/2, L/2, r, -Pi/2, Pi/2, 2*Pi};

f1() = BooleanDifference{ Volume{1}; }{ Volume{2:10}; };
f2() = BooleanIntersection{ Volume{1}; Delete; }{ Volume{2:10}; Delete; };

Physical Volume(0) = {4, 8, 10, 3, 6, 5, 9, 2, 7};
Physical Volume(1) = {11};


Mesh 3;
Coherence Mesh;
