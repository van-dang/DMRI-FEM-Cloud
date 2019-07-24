import os, sys, shutil
from dolfin import *

comm = MPI.comm_world
nprocs = comm.Get_size()
if (nprocs>1):
    print('Support serial computation only!')
    sys.exit()


exists = os.path.isfile('DmriFemLib.py')
isupdate = False
if (exists==False or isupdate==True):
    if isupdate==True:
        os.system("rm DmriFemLib.py")
    print("Load pre-defined functions from GitHub")
    os.system("wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/DmriFemLib.py")
from DmriFemLib import *


"""# Working on the mesh"""

geo_choice = 4
################################################################################
############## Create two-layered disk using mshr in FEniCS ####################
if geo_choice == 1:
    R1, R2 = 5, 10;
    origin = Point(0.,0.)
    circle = Circle(origin, R1, segments=32)
    domain = Circle(origin, R2, segments=32)
    domain.set_subdomain(1, circle)
    mymesh = generate_mesh(domain, 15) # 15 is the resolution
    cmpt_mesh = generate_mesh(circle, 15)
    phase, partion_list, partition_marker = CreatePhaseFunc(mymesh, evengroup, oddgroup, None)

################################################################################
############## Create multilayered domains using gmsh ##########################
if geo_choice == 2:
    mesh_name = "multi_layered_disk"
    # mesh_name = "multi_layered_cylinder"
    # mesh_name = "multi_layered_sphere"
    # mesh_name = "multi_layered_torus"
    
    is_geo_file_exist = os.path.isfile(mesh_name+'.geo')  
    if is_geo_file_exist==False:
        os.system('wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/mesh/'+mesh_name+'.geo')

    # Modify .geo file from 4 layers to 3 layers      
    os.system("sed -i 's/5, 7.5, 10, 13/5, 7.5, 10/g' "+mesh_name+".geo")
      
    # Create mesh from geo file by gmsh
    os.system('gmsh -3 '+mesh_name+'.geo -o '+mesh_name+'.msh')
    
    # Convert .msh to .xml using dolfin-convert
    os.system('dolfin-convert '+mesh_name+'.msh '+mesh_name+'.xml')
       
    mymesh = Mesh(mesh_name+".xml");  

    GetPartitionMarkers(mesh_name+".msh", "pmk_"+mesh_name+".xml")

    partition_marker = MeshFunction("size_t", mymesh, mymesh.topology().dim())

    File("pmk_"+mesh_name+".xml")>>partition_marker

    phase, partion_list = CreatePhaseFunc(mymesh, [], [], partition_marker)    

################################################################################
############## Download the existing mesh and submesh ##########################
if geo_choice == 3:
    is_file_exist = os.path.isfile("multi_layer_torus.xml")  
    if is_file_exist==False:
        os.system('wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/comri/meshes/multi_layer_torus.xml.zip')
        os.system('wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/comri/meshes/multi_layer_torus_compt1.xml.zip')
        os.system('unzip -q multi_layer_torus.xml.zip')
        os.system('unzip -q multi_layer_torus_compt1.xml.zip')

    mymesh = Mesh("multi_layer_torus.xml");  
    cmpt_mesh = Mesh('multi_layer_torus_compt1.xml')
    phase, partion_list, partition_marker = CreatePhaseFunc(mymesh, [], [cmpt_mesh], None)

    
if geo_choice == 4:
    print("Working with volume_box_N_18_7_3_5L_fine.xml")
    is_file_exist = os.path.isfile("volume_box_N_18_7_3_5L_fine.xml")
    if is_file_exist==False:
        os.system('wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/volume_box_N_18_7_3_5L_fine.xml.zip')
        os.system('wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/volume_N_18_7_3_5L_fine.xml.zip')
        os.system('unzip -q volume_box_N_18_7_3_5L_fine.xml.zip')
        os.system('unzip -q volume_N_18_7_3_5L_fine.xml.zip')
    mymesh = Mesh("volume_box_N_18_7_3_5L_fine.xml");
    cmpt_mesh = Mesh('volume_N_18_7_3_5L_fine.xml')
    phase, partion_list, partition_marker = CreatePhaseFunc(mymesh, [], [cmpt_mesh], None)

if geo_choice == 5:
    print("Working with mesh_226cylinders.xml")
    is_mesh_file_exist = os.path.isfile('mesh_226cylinders.xml.zip')
    if is_mesh_file_exist==False:
        os.system('wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/mesh_226cylinders.xml.zip')
	os.system('wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/submesh_226cylinders.xml.zip')
        os.system('rm -rf mesh_226cylinders.xml submesh_226cylinders.xml __MACOSX')
        os.system('unzip mesh_226cylinders.xml.zip')
        os.system('unzip submesh_226cylinders.xml.zip')
    mymesh = Mesh("mesh_226cylinders.xml");
    cmpt_mesh = Mesh('submesh_226cylinders.xml')
    phase, partion_list, partition_marker = CreatePhaseFunc(mymesh, [], [cmpt_mesh], None)
    
################################################################################
############## Save, Plot phase functions and submeshes to verify ##############
print("Save phase function")
File("phase.pvd")<<phase    

print("Partition markers:", partion_list)

V_DG = FunctionSpace(mymesh, 'DG', 0)

D0 = 3e-3
D0_array = [D0, D0]
IC_array = [1, 1]
T2_array = [1e6, 1e6]

# Variable tensor
dofmap_DG = V_DG.dofmap()
d00 = Function(V_DG); d01 = Function(V_DG); d02 = Function(V_DG)
d10 = Function(V_DG); d11 = Function(V_DG); d12 = Function(V_DG)
d20 = Function(V_DG); d21 = Function(V_DG); d22 = Function(V_DG)
T2  = Function(V_DG); disc_ic = Function(V_DG);
        
for cell in cells(mymesh):
      cell_dof = dofmap_DG.cell_dofs(cell.index())
      cmk = partition_marker[cell.index()]
      T2.vector()[cell_dof]      = T2_array[cmk];
      disc_ic.vector()[cell_dof] = IC_array[cmk]; 
      d00.vector()[cell_dof]     = D0_array[cmk]; 
      d11.vector()[cell_dof]     = D0_array[cmk]; 
      d22.vector()[cell_dof]     = D0_array[cmk];

ofile = 'files.h5';

for i in range(0, len(sys.argv)):
      arg = sys.argv[i];
      if arg=='-o':
            ofile = sys.argv[i+1];

filename, file_extension = os.path.splitext(ofile)

ofile = filename+'.h5'

f = HDF5File(mymesh.mpi_comm(), ofile, 'w')
f.write(mymesh, 'mesh');  f.write(T2, 'T2'); f.write(disc_ic, 'ic'); f.write(phase, 'phase');
f.write(d00, 'd00'); f.write(d01, 'd01'); f.write(d02, 'd02')
f.write(d10, 'd10'); f.write(d11, 'd11'); f.write(d12, 'd12')
f.write(d20, 'd20'); f.write(d21, 'd21'); f.write(d22, 'd22')

print("Write to ", ofile)
