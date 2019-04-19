i = 0;
N = 6;
lc = 2.0;
L = 20.0;
xc[]={0,-L,0};
Point(1) = {xc[0], xc[1], xc[2], 	lc};
R[] = {5,7.5,10.0};
ncomp = #R[];
nlayers = 20;
For comp In {0:(ncomp-1)}
  r = R[comp];
  Point(2+i) = {xc[0], r+xc[1], xc[2], 	lc};
  Point(3+i) = {xc[0], xc[1], r+xc[2], 	lc};
  Point(4+i) = {xc[0], xc[1], -r+xc[2], 	lc};
  Point(5+i) = {xc[0], -r+xc[1], xc[2], 	lc};
  Circle(1+i) = {2+i, 1, 4+i};
  Circle(2+i) = {4+i, 1, 5+i};
  Circle(3+i) = {5+i, 1, 3+i};
  Circle(4+i) = {3+i, 1, 2+i};
  Line Loop(1+i) = {1+i, 2+i, 3+i, 4+i};
  i += N;
EndFor
Plane Surface(1) = {1};out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{1}; Layers{nlayers}; };
out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{out[0]}; Layers{nlayers}; };
out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{out[0]}; Layers{nlayers}; };
i = 0;
For r In {1:(ncomp-1)}
  Plane Surface(2+i) = {1+i*N, 1+(i+1)*N};
  out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{2+i}; Layers{nlayers}; };
  out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{out[0]}; Layers{nlayers}; };
  out[] = Extrude { { 0,0,1 }, { 0,0,0 }, 2*Pi/3 } { Surface{out[0]}; Layers{nlayers}; };
  i += 1;
EndFor
Mesh 3;
Coherence Mesh;
