import geompy
import re
import salome
import math
import random
import struct;
import smesh, SMESH, SALOMEDS
import NETGENPlugin
import inspect, os
import libSMESH_Swig

# MultiProjection function : project a set of d-1 mesh on another set of faces
# input : 
#    mesh :        the global mesh
#    source :      a set of faces which lean on a same plane
#    destination : a set of faces which lean on a translated plane 
#                  (the ref plane is the first one)
# -----------------------------------
#def MultiProjection(mesh,source,destination):
#     ls = geompy.SubShapeAllSortedCentres(source     , geompy.ShapeType["EDGE"])
#     ld = geompy.SubShapeAllSortedCentres(destination, geompy.ShapeType["EDGE"])
#     for i in xrange(len(ls)):
#         a = mesh.Projection1D(ld[i])
#         a.SourceEdge(ls[i])

def MultiProjection(probdim,mesh, source,destination):
     if (probdim==2):
	a = mesh.Projection1D(destination)
	a.SourceEdge(source)
     if (probdim==3):
	a = mesh.Projection2D(destination)
	a.SourceFace(source)

# FindFaces function : find in a shape a set of shape leaning on a defined plane
# input : 
#    shape          : the shape
#    origin, normal : definition of the ref plane
#    name           : the name of the set of founded faces
# output:
#    gr : group of the founded faces
# -----------------------------------
def FindFaces(probdim,shape,origin,normal,name):
     if (probdim==2):   
	gr = geompy.CreateGroup(shape, geompy.ShapeType["EDGE"])
	gr.SetName(name)
	geompy.addToStudyInFather(shape, gr, name)
	list = geompy.GetShapesOnPlaneWithLocation(shape, 
		geompy.ShapeType["EDGE"], normal, origin, geompy.GEOM.ST_ON)
     if (probdim==3):   
	gr = geompy.CreateGroup(shape, geompy.ShapeType["FACE"])
	gr.SetName(name)
	geompy.addToStudyInFather(shape, gr, name)
	list = geompy.GetShapesOnPlaneWithLocation(shape, 
		geompy.ShapeType["FACE"], normal, origin, geompy.GEOM.ST_ON)
     geompy.UnionList(gr, list)       
     return gr

directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/' # script directory
xstart = -5;
ystart = -5;
zstart = -5;
xend = 5;
yend = 5;
zend = 5;

Box = geompy.MakeBox(xstart,ystart,zstart,xend,yend,zend);
Box_file = directory+"Box.brep"
TrimSizePlane=1000.0;
origin = geompy.MakeVertex(0, 0, 0)
vx = geompy.MakeVectorDXDYDZ(1, 0, 0)
vy = geompy.MakeVectorDXDYDZ(0, 1, 0)
vz = geompy.MakeVectorDXDYDZ(0, 0, 1)
planeX = geompy.MakePlane(origin,vx,TrimSizePlane)
planeY = geompy.MakePlane(origin,vy,TrimSizePlane)
planeZ = geompy.MakePlane(origin,vz,TrimSizePlane)
planes = [planeX,planeY,planeZ]
Box = geompy.MakePartition([Box],planes);
geompy.addToStudy(Box,"Box")

center=[[0,0,0], [-5,-5,-5], [-5,5,-5], [5,-5,-5], [5,5,-5], [-5,-5,5], [-5,5,5], [5,-5,5], [5,5,5]]

Rmean=4.0;


ArrayOfSpheres1=[];

for i in xrange(0,len(center)):
  x0=center[i][0]; y0=center[i][1]; z0=center[i][2]; 
  ArrayOfSpheres1.append(geompy.MakeSphere(x0,y0,z0,Rmean));


GroupOfSphere1 = geompy.MakeCompound( ArrayOfSpheres1);
GroupOfSphere1 = geompy.MakeCommon(GroupOfSphere1, Box);

geompy.addToStudy(GroupOfSphere1,'GroupOfSphere1');


Group1 = GroupOfSphere1;


geompy.addToStudy(Group1 ,'Group1');


FinalPartition = geompy.MakePartition([Box,Group1],[], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
geompy.addToStudy( FinalPartition, 'FinalPartition' )


Gp1 = geompy.CreateGroup(FinalPartition, geompy.ShapeType["SOLID"])
geompy.addToStudyInFather(FinalPartition, Gp1, 'Gp1' )
subGp1 = geompy.GetInPlace(FinalPartition,Group1,True)
geompy.UnionList(Gp1, [subGp1])


vx = geompy.MakeVectorDXDYDZ(1, 0, 0)
vy = geompy.MakeVectorDXDYDZ(0, 1, 0)
vz = geompy.MakeVectorDXDYDZ(0, 0, 1)

le = FindFaces(3,FinalPartition,geompy.MakeVertex(0., ystart, 0.),vy,"left")
ri = FindFaces(3,FinalPartition,geompy.MakeVertex(0.,  yend, 0.),vy,"right")
fr = FindFaces(3,FinalPartition,geompy.MakeVertex( xend, 0., 0.),vx,"front")
be = FindFaces(3,FinalPartition,geompy.MakeVertex(xstart, 0., 0.),vx,"behind")
bo = FindFaces(3,FinalPartition,geompy.MakeVertex(0., 0.,  zend),vz,"bottom")
to = FindFaces(3,FinalPartition,geompy.MakeVertex(0., 0., zstart),vz,"top")

FinalMesh = smesh.Mesh(FinalPartition,"FinalPartition")
algo = FinalMesh.Tetrahedron(smesh.NETGEN)
algo1d = FinalMesh.Segment()
algo1d.LocalLength(0.5)

algo2d = FinalMesh.Triangle(smesh.NETGEN_2D)


MultiProjection(3,FinalMesh, le, ri)
MultiProjection(3,FinalMesh, fr, be)
MultiProjection(3,FinalMesh, bo, to)

FinalMesh.Compute()

tol = 1e-4;
GroupsOfNodes = FinalMesh.FindCoincidentNodes(tol);
FinalMesh.MergeNodes(GroupsOfNodes);

import SMESH
GroupMesh=FinalMesh.GroupOnGeom(Gp1,"1",SMESH.VOLUME)

# Step 7 : Save mesh
# -----------------------------------
FinalMesh.ExportUNV(directory+"SpheresInBox.unv")
FinalMesh.ExportUNV( directory+"SpheresInBox_cmpt1.unv", meshPart=GroupMesh )
