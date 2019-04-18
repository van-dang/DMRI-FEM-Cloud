lc = 1.0;
i = 0;
N = 30;
R[] = {5, 7.5, 10};
ncomp = #R[];
Printf("Number compartments: %g", ncomp);
For comp In {0:(ncomp-1)}
    r = R[comp];
    Printf("Working on compartment %g", comp);
    Point(1+i) = {0.0,0.0,0.0,lc};
    Point(2+i) = {r,0.0,0.0,lc};
    Point(3+i) = {0,r,0.0,lc};
    Circle(1+i) = {2+i,1+i,3+i};
    Point(4+i) = {-r,0,0.0,lc};
    Point(5+i) = {0,-r,0.0,lc};
    Circle(2+i) = {3+i,1+i,4+i};
    Circle(3+i) = {4+i,1+i,5+i};
    Circle(4+i) = {5+i,1+i,2+i};
    Point(6+i) = {0,0,-r,lc};
    Point(7+i) = {0,0,r,lc};
    Circle(5+i) = {3+i,1+i,6+i};
    Circle(6+i) = {6+i,1+i,5+i};
    Circle(7+i) = {5+i,1+i,7+i};
    Circle(8+i) = {7+i,1+i,3+i};
    Circle(9+i) = {2+i,1+i,7+i};
    Circle(10+i) = {7+i,1+i,4+i};
    Circle(11+i) = {4+i,1+i,6+i};
    Circle(12+i) = {6+i,1+i,2+i};
    Line Loop(13+i) = {2+i,8+i,-(10+i)};
    Ruled Surface(14+i) = {13+i};
    Line Loop(15+i) = {10+i,3+i,7+i};
    Ruled Surface(16+i) = {15+i};
    Line Loop(17+i) = {-(8+i),-(9+i),1+i};
    Ruled Surface(18+i) = {17+i};
    Line Loop(19+i) = {-(11+i),-(2+i),5+i};
    Ruled Surface(20+i) = {19+i};
    Line Loop(21+i) = {-(5+i),-(12+i),-(1+i)};
    Ruled Surface(22+i) = {21+i};
    Line Loop(23+i) = {-(3+i),11+i,6+i};
    Ruled Surface(24+i) = {23+i};
    Line Loop(25+i) = {-(7+i),4+i,9+i};
    Ruled Surface(26+i) = {25+i};
    Line Loop(27+i) = {-(4+i),12+i,-(6+i)};
    Ruled Surface(28+i) = {27+i};
    Surface Loop(29+i) = {28+i,26+i,16+i,14+i,20+i,24+i,22+i,18+i};
    i = i + N;
EndFor

Volume(1) = {29};
j = 0;
For comp In {1:(ncomp-1)}
    Volume(j+2) = {29+j*N, 29+(j+1)*N};
    j = j + 1;
EndFor

Mesh 3;
Coherence Mesh; // Save "spherical_layers.msh";