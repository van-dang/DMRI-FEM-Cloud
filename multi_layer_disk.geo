lc = 1.0;
R[] = {5, 7.5, 10, 13};
ncomp = #R[];
Printf("Number compartments %g", ncomp);
jump=0;
j=0;
For comp In {0:ncomp-1}
	Printf("compartment %g, radius %g", comp, R[comp]);
	r=R[comp];
	Point(1+jump) = {-r, 0, 0, lc};
	Point(2+jump) = {r, 0, 0, lc};
	Point(3+jump) = {0, 0, 0, lc};
	Circle(1+jump) = {1+jump,3+jump,2+jump};
	Circle(2+jump) = {2+jump,3+jump,1+jump};
	If (jump==0) 
	   	Line Loop(10) = {1, 2};
	Else
		Line Loop(10+jump) = {1+jump, 2+jump, 1+jump-3, 2+jump-3};
	EndIf
	Plane Surface(j++) = {10+jump};
	jump=jump+3;
EndFor
Mesh 2;
Coherence Mesh;
Printf("Completed the task!");
