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


# Copyright Dang Van NGUYEN 2015

from MEDLoader import *
from MEDCoupling import *
import os, inspect

# MultiProjection function : project a set of d-1 mesh on another set of 
# faces
# input : 
#    mesh :        the global mesh
#    source :      a set of faces which lean on a same plane
#    destination : a set of faces which lean on a translated plane 
#                  (the ref plane is the first one)
# -----------------------------------
def MultiProjection(probdim,mesh,source,destination,is_cell_cut_bc=False):
  if (is_cell_cut_bc==True):
    # this method is better when cells cut the external boundaries
    if (probdim==2):
      ls = geompy.SubShapeAllSortedCentres(source     , geompy.ShapeType["EDGE"])
      ld = geompy.SubShapeAllSortedCentres(destination, geompy.ShapeType["EDGE"])
      for i in xrange(len(ls)):
	a = mesh.Projection1D(ld[i])
	a.SourceEdge(ls[i])
    if (probdim==3):
      ls = geompy.SubShapeAllSortedCentres(source     , geompy.ShapeType["FACE"])
      ld = geompy.SubShapeAllSortedCentres(destination, geompy.ShapeType["FACE"])
      for i in xrange(len(ls)):
        a = mesh.Projection2D(ld[i])
	a.SourceFace(ls[i])
  else:
    if (probdim==2):
	a = mesh.Projection1D(destination)
	a.SourceEdge(source)
    if (probdim==3):
	a = mesh.Projection2D(destination)
	a.SourceFace(source)

# FindFaces function : find in a shape a set of shape leaning on a defined
# plane
# input : 
#    shape          : the shape
#    origin, normal : definition of the ref plane
#    name           : the name of the set of found faces
# output:
#    gr : group of the found faces
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

def unv2xml(directory, unvfilename, is_use_group_name=True,is_regroup=False, markername=[]):
	# directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/' # script directory
	unvFile = directory + unvfilename;
	unvmesh = smesh.CreateMeshesFromUNV(unvFile);

	filewithoutextention = os.path.splitext(os.path.basename(unvFile))[0];

	medFile = directory+filewithoutextention+'.med';
	unvmesh.ExportMED(medFile)

	# medmesh = unvmesh.GetMEDMesh()
	medmesh = MEDFileMesh.New(medFile)

	print "Num of vertices: ",medmesh.getNumberOfNodes()
	dim = medmesh.getMeshDimension()
	print "Mesh dimension: ",dim

	numvertices = medmesh.getNumberOfNodes();
	if dim==2:
	  numcells = unvmesh.NbTriangles()
	if dim==3:
	  numcells = unvmesh.NbVolumes()
	print "Num of cells: ", numcells;

	# Step 2 : Retrieve the lists of Ids of the cells
	# -----------------------------------

	groupnames=medmesh.getGroupsNames();
	num_groups = len(groupnames)
	print 'num groups: ', num_groups
	if (num_groups>0):
	  if is_use_group_name == True:
	     markerarray = [0]*numcells
	  else:
	     markerarray = [num_groups]*numcells
	     # markername = list(xrange(num_groups));
	  for i in xrange(0,num_groups):
		groupids = medmesh.getGroupArr(0,groupnames[i]).getValues()
		if len(groupids)==0:
		  raise ValueError("Error: Groupe "+ groupnames[i]+" is empty!")
		for k in xrange(0,len(groupids)):
			# data = "%d\n" %groupids[k]; filename.write(data);
			if is_use_group_name == True:
			  name = int(groupnames[i]);
			  if is_regroup==True:
			    name = markername[int(groupnames[i])];
			  markerarray[groupids[k]] = name;
			else:
			  name = markername[i]
			  markerarray[groupids[k]] = name; 
			if k==0:
			  print 'Group', groupnames[i],' is marked with marker', name

	  filename = open(directory+'theta_'+filewithoutextention+'.xml','w');
	  filename.write('<?xml version="1.0"?>\n')
	  filename.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
	  filename.write('  <mesh_function>\n')
	  filename.write('    <mesh_value_collection type="uint" dim="'+str(dim)+'" size="'+str(numcells)+'">\n');
	  for i in xrange(0,len(markerarray)):
	    filename.write('      <value cell_index="'+str(i)+'" local_entity="0" value="'+str(markerarray[i])+'" />\n');
	  filename.write('    </mesh_value_collection>\n');
	  filename.write('  </mesh_function>\n');
	  filename.write('</dolfin>\n');
	  filename.close();

	os.system('rm '+medFile);

	# export .xml mesh
	filename = open(directory+filewithoutextention+'.xml','w');
	filename.write('<?xml version="1.0"?>\n')
	filename.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
	if dim==2:
		filename.write('  <mesh celltype="triangle" dim="2">\n');
	if dim==3:
		filename.write('  <mesh celltype="tetrahedron" dim="3">\n');
	# save mesh vertices
	filename.write('   <vertices size="'+str(numvertices)+'">\n');

	for node in xrange(1,numvertices+1):
	  xyz = unvmesh.GetNodeXYZ( node );
	  if (dim==2):		
	   	filename.write('	<vertex index="'+str(node-1)+'" x="%.16e" y="%.16e" />\n' % (xyz[0],xyz[1]));	
	  if (dim==3):		
	   	filename.write('	<vertex index="'+str(node-1)+'" x="%.16e" y="%.16e" z="%.16e" />\n' % (xyz[0],xyz[1],xyz[2]));	
	filename.write('   </vertices>\n');

	# save mesh cells
	filename.write('   <cells size="'+str(numcells)+'">\n');
	if dim==2:
		vol = unvmesh.GetElementsByType(SMESH.FACE) 
	if dim==3:
		vol = unvmesh.GetElementsByType(SMESH.VOLUME) 
	# print vol
	cellindex = 0;
	for ind in vol:
	  nodesid = unvmesh.GetElemNodes(ind)
	  if (dim==2):
		filename.write('	<triangle index="%d" v0="%d" v1="%d" v2="%d" />\n' % (cellindex,nodesid[0]-1,nodesid[1]-1,nodesid[2]-1) );
	  if (dim==3):
		filename.write('	<tetrahedron index="'+str(cellindex)+'" v0="'+str(nodesid[0]-1)+'" v1="'+str(nodesid[2]-1)+'" v2="'+str(nodesid[1]-1)+'" v3="'+str(nodesid[3]-1)+'" />\n');
	  cellindex = cellindex + 1;

	filename.write('   </cells>\n');	

	filename.write('  </mesh>\n');
	filename.write('</dolfin>\n');
	filename.close();
	if (num_groups>0):
	  print "Program is terminal normally. Two files were created:";
	  print 'unvmesh:   '+ directory+filewithoutextention+'.xml';
	  print 'makerfile: '+directory+'theta_'+filewithoutextention+'.xml';
	else:
	  print "Program is terminal normally. One file was created:";
	  print 'unvmesh:   '+ directory+filewithoutextention+'.xml';


def checkProjection(gr, mesh_translated, tol=1e-7):
    name = gr.GetName() + "_" + mesh_translated.GetName().split("_")[0]
    mesh_source = smesh.CopyMesh(gr, gr.GetName())
    mesh_check = smesh.Concatenate([mesh_source.GetMesh(), mesh_translated.GetMesh()], 0)
    mesh_check.SetName(name)
    ll_coincident_nodes = mesh_check.FindCoincidentNodes(tol)
    coincident_nodes = [item for sublist in ll_coincident_nodes for item in sublist]
    mesh_check.MakeGroupByIds("coincident_nodes", smesh.NODE, coincident_nodes)
    mesh_nodes = mesh_check.GetNodesId()
    if len(ll_coincident_nodes) != mesh_translated.NbNodes():
        print "Projection failed for %s"%name
        non_coincident_nodes = list(set(mesh_nodes) - set(coincident_nodes))
        mesh_check.MakeGroupByIds("non_coincident_nodes", smesh.NODE, non_coincident_nodes)


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


MultiProjection(3,FinalMesh, le, ri, True)
MultiProjection(3,FinalMesh, fr, be, True)
MultiProjection(3,FinalMesh, bo, to, True)

FinalMesh.Compute()

tol = 1e-4;
GroupsOfNodes = FinalMesh.FindCoincidentNodes(tol);
FinalMesh.MergeNodes(GroupsOfNodes);

# Test the peridic boundaries 
gr_le = FinalMesh.Group(le) 
gr_ri = FinalMesh.Group(ri) 
gr_be = FinalMesh.Group(be) 
gr_fr = FinalMesh.Group(fr) 
gr_bo = FinalMesh.Group(bo) 
gr_to = FinalMesh.Group(to)

le_translated = FinalMesh.TranslateObjectMakeMesh( gr_le, smesh.DirStruct(smesh.PointStruct ( 0, yend -ystart , 0 )), 0, 'le_translated')
be_translated = FinalMesh.TranslateObjectMakeMesh( gr_be, smesh.DirStruct(smesh.PointStruct ( -xstart+xend, 0, 0 )), 0, 'be_translated')
bo_translated = FinalMesh.TranslateObjectMakeMesh( gr_bo, smesh.DirStruct(smesh.PointStruct ( 0, 0, -zend+zstart )), 0, 'bo_translated')

checkProjection(gr_ri, le_translated) 
checkProjection(gr_fr, be_translated) 
checkProjection(gr_to, bo_translated) 

import SMESH
GroupMesh=FinalMesh.GroupOnGeom(Gp1,"1",SMESH.VOLUME)

# Step 7 : Save mesh
# -----------------------------------
FinalMesh.ExportUNV(directory+"SpheresInBox.unv")
FinalMesh.ExportUNV( directory+"SpheresInBox_cmpt1.unv", meshPart=GroupMesh )
