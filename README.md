# ABOUT

DMRI-FEM-Cloud is an effort to bring Cloud Computing to computational diffusion MRI. It supports

(1) Single domains, Multilayered structures, Manifolds
(2) Membrane permeability for internal interfaces and Artificial permeability at the external interfaces
(3) Pure homogeneous Neumann BCs and Pseudo-periodic BCs
(4) Web interface Google Colab / Python notebooks and TUI.
(5) Serial and parallel executions

Copyright (C) 2019 Van-Dang Nguyen

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

# Simulations with Python notebooks

We designed some benchmark problems that illustrate how the diffusion MRI simulations can be performed on Google Colab. Users can open the files in Github and execute it directly with Google Colab. The list of benchmarks problems are as follows

* [Arbitrary diffusion pulses](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/ArbitraryTimeSequence.ipynb) (PGSE, OGSE, Double-PGSE, Double-OGSE, Trapezoidal PGSE,...)

* [Extracellular Space](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/ExtracellularSpace.ipynb)

* [Realistic neurons](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/RealNeurons.ipynb)

* [Manifolds](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/Manifolds.ipynb)

* [Discontinuous Initial Conditions](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/DiscontinuousInitialCondition.ipynb)

* [Multilayered Structures](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/MultilayeredStructures.ipynb)

* [Multilayered Disk with Variable Permeability](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/MultilayeredDiskVariablePermeability.ipynb)

* [Pseudo-Periodicity](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/PeriodicDomains.ipynb)

* [T2-relaxation](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/T2_Relaxation.ipynb)

* [Explicit Implementation](https://colab.research.google.com/github/van-dang/DMRI-FEM-Cloud/blob/master/ExplicitImplementation.ipynb)

These notebooks can be connected to either a hosted runtime provided by Google Cloud or a local runtime. The hosted runtime allows us to access free resources for up to 12 hours at a time. For longer excecutions, it is more convenient to connect to the local runtimes. The instructions are available at [this link](https://github.com/van-dang/DMRI-FEM-Cloud/blob/master/LocalColab.md).

# Simulations with Singularity images on Google Cloud

## Create a VM instance

Access https://cloud.google.com

Go to console

Nivigation menu / Compute Engine / VM instances / Create Instances / SSH connect

```bash
sudo apt-get install singularity-container unzip
```

## With FEniCS
### Create a FEniCS Image in writable mode

```bash
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/singularity_images/Singularity_recipe_FEniCS_DMRI
sudo singularity build -w writable_fenics_dmri.simg Singularity_recipe_FEniCS_DMRI
```

### Test if mpi4py works correctly
```bash
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/test_mpi4py.py
mpirun -n 3 singularity exec -B $PWD writable_fenics_dmri.simg python3 test_mpi4py.py
```
The results would be
```bash
My rank is  1
My rank is  2
My rank is  0
```

### Copy Python solvers to the VM instance
```bash
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/PreprocessingOneCompt.py
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/PreprocessingMultiCompt.py
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/GCloudDmriSolver.py
```

### For multi-compartment domains
```bash
singularity exec -B $PWD writable_fenics_dmri.simg python3 PreprocessingMultiCompt.py -o multcompt_files.h5
mpirun -n 8 singularity exec -B $PWD writable_fenics_dmri.simg python3 GCloudDmriSolver.py -f multcompt_files.h5 -M 1 -b 1000 -p 1e-5 -d 10600 -D 43100 -k 200 -gdir 0 1 0
 ```
### For single-compartment domains
```bash
singularity exec -B $PWD writable_fenics_dmri.simg python3 PreprocessingOneCompt.py -o onecompt_files.h5
mpirun -n 8 singularity exec -B $PWD writable_fenics_dmri.simg python3 GCloudDmriSolver.py -f onecompt_files.h5 -M 0 -b 1000 -d 10600 -D 43100 -k 200 -K 3e-3 -gdir 1 0 0 
```
## With FEniCS-HPC
### Create a FEniCS Image in writable mode

```bash
wget https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/singularity_images/Singularity_recipe_FEniCS_HPC_DMRI
sudo singularity build -w writable_fenics_hpc_dmri.simg Singularity_recipe_FEniCS_HPC_DMRI
```

### Test if mpi works correctly
```bash
wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
singularity exec -B $PWD writable_fenics_hpc_dmri.simg mpicc mpi_hello_world.c -o mpi_hello_world
mpirun -n 3 singularity exec -B $PWD writable_fenics_hpc_dmri.simg  mpi_hello_world
```
The results would be
```bash
Hello world from processor dmri, rank 0 out of 3 processors
Hello world from processor dmri, rank 1 out of 3 processors
Hello world from processor dmri, rank 2 out of 3 processors
```

### Download the solvers
```bash
wget https://github.com/van-dang/DMRI-FEM-Cloud/archive/fenics-hpc-solvers.zip
unzip fenics-hpc-solvers.zip
```

### For single-compartment domains

```bash
# Compile the form files
cd DMRI-FEM-Cloud-fenics-hpc-solvers/one-comp/ufc
singularity exec -B $PWD ../../../writable_fenics_hpc_dmri.simg make -j 8
cd ../

# Compile main.cpp
singularity exec -B $PWD ../../writable_fenics_hpc_dmri.simg make clean
singularity exec -B $PWD ../../writable_fenics_hpc_dmri.simg make -j 8

# Create a working directory
mkdir test_04b_pyramidal7aACC
cd test_04b_pyramidal7aACC

# Copy the executable demo to the working directory
cp ../demo .

# Download the mesh
wget https://github.com/van-dang/RealNeuronMeshes/raw/master/volume_meshes/pyramidals/04b_pyramidal7aACC.msh.zip
unzip 04b_pyramidal7aACC.msh.zip
# Convert .msh to .xml
wget https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.py
python dolfin-convert.py 04b_pyramidal7aACC.msh 04b_pyramidal7aACC.xml

# Execute the demo
mpirun -n 8 singularity exec -B $PWD ../../../writable_fenics_hpc_dmri.simg ./demo -m 04b_pyramidal7aACC.xml -b 1000 -d 10600 -D 43100 -k 200 -K 3e-3 -v 1 0 0  > my_output_file
```

### For two-compartment domains
```bash
# Compile the form files
cd DMRI-FEM-Cloud-fenics-hpc-solvers/two-comp/ufc
singularity exec -B $PWD ../../../writable_fenics_hpc_dmri.simg make -j 8
cd ../

# Compile main.cpp
singularity exec -B $PWD ../../writable_fenics_hpc_dmri.simg make clean
singularity exec -B $PWD ../../writable_fenics_hpc_dmri.simg make -j 8

# Create a working directory
mkdir test_neuron_N_18_7_3_5L
cd test_neuron_N_18_7_3_5L

# Copy the executable demo to the working directory
cp ../demo .

# Download the existing meshes
wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/volume_box_N_18_7_3_5L_fine.xml.zip
wget --quiet https://github.com/van-dang/DMRI-FEM-Cloud/raw/mesh/volume_N_18_7_3_5L_fine.xml.zip
unzip -q volume_box_N_18_7_3_5L_fine.xml.zip
unzip -q volume_N_18_7_3_5L_fine.xml.zip

# Execute the demo
mpirun -8 singularity exec -B $PWD ../../../writable_fenics_hpc_dmri.simg  ./demo -m volume_box_N_18_7_3_5L_fine.xml -c volume_N_18_7_3_5L_fine.xml -b 1000 -p 1e-5 -d 10600 -D 43100 -k 200 -v 1 0 0 
```

