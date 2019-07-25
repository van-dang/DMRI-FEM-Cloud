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

geo_choice = 2

if geo_choice == 0:
      mesh_file='2E_ExtraCellular_group_10um_vol'
      file_dir='https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/mesh/'+mesh_file+'.msh.zip'
      zip_exists = os.path.isfile(mesh_file+".msh.zip")
      mesh_file_exists = os.path.isfile(mesh_file)
      if (zip_exists==False):
            os.system("wget "+file_dir)
      if (mesh_file_exists==False):
            os.system("unzip -q "+mesh_file+".msh.zip")
      os.system("dolfin-convert "+mesh_file+".msh "+mesh_file+".xml")
                   
if geo_choice == 1:
      mesh_file = "fru_M_100383_1D.xml"
      zip_exists = os.path.isfile(mesh_file+".zip")
      mesh_file_exists = os.path.isfile(mesh_file)
      if (zip_exists==False):
            os.system("wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/mesh/"+mesh_file+".zip")
      if (mesh_file_exists==False):
            os.system("unzip -q "+mesh_file+".zip")
      mymesh = Mesh(mesh_file);

if geo_choice == 2:
      spindle_list=['03a_spindle2aFI','03a_spindle6aFI','03b_spindle4aACC','03b_spindle5aACC',
                    '03b_spindle6aACC','03b_spindle7aACC','04b_spindle3aFI','05b_spindle5aFI',
                    '06b_spindle8aACC','07b_spindle9aACC','08a_spindle13aACC','09o_spindle7aFI',
                    '09o_spindle8aFI','10a_spindle18aACC','12a_spindle19aACC','12o_spindle9aFI',
                    '13o_spindle10aFI','15o_spindle12aFI','16o_spindle13aFI','19o_spindle14aFI',
                    '21o_spindle15aFI','23o_spindle16aFI','25o_spindle17aFI','26o_spindle18aFI',
                    '27o_spindle19aFI','28o_spindle20aFI','28o_spindle21aFI','29o_spindle22aFI',
                    '30o_spindle23aFI',
      ];
      pyramidal_list=['02a_pyramidal2aFI','02b_pyramidal1aACC','02b_pyramidal1aFI','03a_pyramidal9aFI',
                      '03b_pyramidal2aACC','03b_pyramidal3aACC','03b_pyramidal3aFI','03b_pyramidal4aFI',
                      '03b_pyramidal9aFI','04a_pyramidal4aACC','04a_pyramidal5aACC','04b_pyramidal5aFI',
                      '04b_pyramidal6aACC','04b_pyramidal6aFI','04b_pyramidal7aACC','05a_pyramidal10aACC',
                      '05a_pyramidal8aACC','05b_pyramidal7aFI','05b_pyramidal8aFI','05b_pyramidal9aACC',
                      '06a_pyramidal11aACC','06b_pyramidal10aFI','06b_pyramidal12aACC','07a_pyramidal13aACC',
                      '07b_pyramidal14aACC','08o_pyramidal11aFI','10a_pyramidal15aACC','11a_pyramidal16aACC',
                      '11o_pyramidal12aFI','17o_pyramidal13aFI','18o_pyramidal14aFI','20o_pyramidal15aFI',
                      '22o_pyramidal16aFI','24o_pyramidal17aFI','25o_pyramidal18aFI','31o_pyramidal19aFI',
      ];

      neuron_id = 14;
      neuron_type = 'pyramidals';

      if neuron_type == 'spindles':
          neuron_list = spindle_list;
      if neuron_type == 'pyramidals':
          neuron_list = pyramidal_list;

      mesh_file = neuron_list[neuron_id];
      
      neuron_dir='https://github.com/van-dang/RealNeuronMeshes/raw/master/volume_meshes/'+neuron_type+'/'+mesh_file+'.msh.zip'

      zip_exists = os.path.isfile(mesh_file+".msh.zip")

      if (zip_exists==False):
              os.system("wget -q "+neuron_dir)
      os.system("unzip -q "+mesh_file+".msh.zip")
      os.system("dolfin-convert "+mesh_file+".msh "+mesh_file+".xml")

print('Mesh file: ', mesh_file)

mesh = Mesh(mesh_file+".xml")

V_DG = FunctionSpace(mesh, 'DG', 0)

D0_array = [3e-3]
IC_array = [1]
T2_array = [1e6]

# Variable tensor
dofmap_DG = V_DG.dofmap()
d00 = Function(V_DG); d01 = Function(V_DG); d02 = Function(V_DG)
d10 = Function(V_DG); d11 = Function(V_DG); d12 = Function(V_DG)
d20 = Function(V_DG); d21 = Function(V_DG); d22 = Function(V_DG)
T2 = Function(V_DG); disc_ic = Function(V_DG);
        

print('Setting parameters to %d cells'%(mesh.num_cells()))

T2.vector()[:]      = T2_array[0];
disc_ic.vector()[:] = IC_array[0];
d00.vector()[:]     = D0_array[0];
d11.vector()[:]     = D0_array[0];
d22.vector()[:]     = D0_array[0];
                             
'''
for cell in cells(mesh):
      cell_dof = dofmap_DG.cell_dofs(cell.index())
      T2.vector()[cell_dof]      = T2_array[0];
      disc_ic.vector()[cell_dof] = IC_array[0];
      d00.vector()[cell_dof]     = D0_array[0];
      d11.vector()[cell_dof]     = D0_array[0];
      d22.vector()[cell_dof]     = D0_array[0];
'''

ofile = 'files.h5';

for i in range(0, len(sys.argv)):
      arg = sys.argv[i];
      if arg=='-o':
            ofile = sys.argv[i+1];

filename, file_extension = os.path.splitext(ofile)


ofile = filename+'.h5'

print("Write to ",ofile)
f = HDF5File(mesh.mpi_comm(), ofile, 'w')
f.write(mesh, 'mesh')
f.write(T2, 'T2');  f.write(disc_ic, 'ic'); 
f.write(d00, 'd00'); f.write(d01, 'd01'); f.write(d02, 'd02')
f.write(d10, 'd10'); f.write(d11, 'd11'); f.write(d12, 'd12')
f.write(d20, 'd20'); f.write(d21, 'd21'); f.write(d22, 'd22')

print("Done")
