# Simulations with Python notebooks

The scope of usage: 
(1) Single domains, Multilayered structures, manifolds
(2) Membrane permeability for internal interfaces
    Artificial permeability at the external interfaces
(3) pure homogeneous Neumann BCs, (4) pseudo-periodic BCs

Copyright (C) 2019 Van-Dang Nguyen

We designed some examples that illustrate how the diffusion MRI simulations can be performed on Google Colab. Users can open the files in Github and execute it directly with Google Colab.

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/ArbitraryTimeSequence.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/ExtracellularSpace.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/RealNeurons.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/Manifolds.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/DiscontinuousInitialCondition.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/MultilayeredStructures.ipynb

https://colab.research.google.com/github/van-dang/MRI-Cloud/blob/master/T2_Relaxation.ipynb

However, if there is any problem with the Github visualization, users can follow two steps:

1. Download the files into local computer: Users can clone the whole repository using 'git clone git@github.com:van-dang/MRI-Cloud.git'. However, it is possible to download the files individually in the raw mode for instance:
'wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/T2_Relaxation.ipynb'

2. Open a new Colab notebook by clicking on https://colab.research.google.com and upload the downloaded notebook 'File/Upload notebook ...'


# Simulations on Google Cloud

## Create a VM instance

Access https://cloud.google.com

Go to console

Nivigation menu / Compute Engine / VM instances / Create Instances / SSH connect

```bash
sudo apt-get install singularity-container unzip
```

## Working with FEniCS Singularity Image
### Create a FEniCS Image in writable mode

```bash
sudo singularity build --writable writable_fenics_stable.simg docker://fenicsproject/stable
```
### Or download the exisiting FEniCS Image and change it to writable mode
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics_stable.simg
sudo singularity build --writable writable_fenics_stable.simg fenics_stable.simg
```

### Install some dependencies
```bash
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get update
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get install zip unzip gmsh
```
### Test if mpi4py works correctly
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/test_mpi4py.py
mpirun -n 3 singularity exec -B $PWD writable_fenics_stable.simg python3 test_mpi4py.py
```
The results would be
```bash
My rank is  1
My rank is  2
My rank is  0
```

### Copy Python solvers to the VM instance
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/PreprocessingOneCompt.py
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/PreprocessingMultiCompt.py
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/GCloudDmriSolver.py
```

### Execute the code with Singularity Image

##### For multi-compartment domains
```bash
singularity exec -B $PWD writable_fenics_stable.simg python3 PreprocessingMultiCompt.py -o multcompt_files.h5
mpirun -n 8 singularity exec -B $PWD writable_fenics_stable.simg python3 GCloudDmriSolver.py -f multcompt_files.h5 -M 1 -b 1000 -p 1e-5 -d 10600 -D 43100 -k 200 -gdir 0 1 0
 ```
##### For single-compartment domains
```bash
singularity exec -B $PWD writable_fenics_stable.simg python3 PreprocessingOneCompt.py -o onecompt_files.h5
mpirun -n 8 singularity exec -B $PWD writable_fenics_stable.simg python3 GCloudDmriSolver.py -f onecompt_files.h5 -M 0 -b 1000 -d 10600 -D 43100 -k 200 -K 3e-3 -gdir 1 0 0 
```
## Working with FEniCS-HPC Singularity Image
### Download existing images
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
```

### Change to writable mode
```bash
sudo singularity build --writable writable_fenics-hpc-dmri.simg fenics-hpc-dmri.simg
```

### Install packages to the existing image
```bash
sudo singularity exec --writable writable_fenics-hpc-dmri.simg apt-get install zip unzip gmsh
```

### Test if mpi works correctly
```bash
wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
singularity exec -B $PWD writable_fenics-hpc-dmri.simg mpicc mpi_hello_world.c -o mpi_hello_world
singularity exec -B $PWD writable_fenics-hpc-dmri.simg mpirun -n 3  mpi_hello_world
```
The results would be
```bash
Hello world from processor dmri, rank 0 out of 3 processors
Hello world from processor dmri, rank 1 out of 3 processors
Hello world from processor dmri, rank 2 out of 3 processors
```

### Download the solvers
```bash
wget https://github.com/van-dang/MRI-Cloud/archive/fenics-hpc-solvers.zip
unzip fenics-hpc-solvers.zip
cd MRI-Cloud-fenics-hpc-solvers/one-comp/
```
### Compile
```bash
singularity exec -B $PWD ../../writable_fenics-hpc-dmri.simg make -j 8
```
### Download the existing meshes
```bash
https://github.com/van-dang/RealNeuronMeshes/raw/master/volume_meshes/pyramidals/04b_pyramidal7aACC.msh.zip
unzip 04b_pyramidal7aACC.msh.zip
```

### Convert .gmsh to .xml
```bash
wget https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.py
python dolfin-convert.py 04b_pyramidal7aACC.msh 04b_pyramidal7aACC.xml
```

### Execute the demo
```bash
singularity exec -B $PWD ../../fenics-hpc-dmri.simg mpirun -n 8 ./demo -m 04b_pyramidal7aACC.xml -b 1000 -d 10600 -D 43100 -k 200 -K 3e-3 -v 1 0 0  > my_output_file
```
