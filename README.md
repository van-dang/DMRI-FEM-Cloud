# Simulations with Python notebooks

The scope of usage: 
(1) Single domains, Multilayered structures, manifolds
(2) Membrane permeability for internal interfaces
    Artificial permeability at the external interfaces
(3) pure homogeneous Neumann BCs, (4) pseudo-periodic BCs

Copyright (C) 2019 Van-Dang Nguyen

We designed some examples that illustrate how the diffusion MRI simulations can be performed on Google Colab. Users can open the files in Github and execute it directly with Google Colab.

https://github.com/van-dang/MRI-Cloud/blob/master/ArbitraryTimeSequence.ipynb

https://github.com/van-dang/MRI-Cloud/blob/master/ExtracellularSpace.ipynb

https://github.com/van-dang/MRI-Cloud/blob/master/RealNeurons.ipynb

https://github.com/van-dang/MRI-Cloud/blob/master/DiscontinuousInitialCondition.ipynb

https://github.com/van-dang/MRI-Cloud/blob/master/MultilayeredStructures.ipynb

https://github.com/van-dang/MRI-Cloud/blob/master/T2_Relaxation.ipynb

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
sudo apt-get update 
sudo apt-get install singularity-container
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
### Test mpi4py the image
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/test_mpi4py.py
mpirun -n 3 singularity exec -B $PWD writable_fenics_stable.simg python3 test.py
```
The results would be
```bash
My rank is  1
My rank is  2
My rank is  0
```

### Install some dependencies
```bash
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get update
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get install zip unzip gmsh
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
singularity exec -B $PWD writable_fenics_stable.simg mpirun -n 6 python3 GCloudDmriSolver.py -f multcompt_files.h5 -M 1 -b 1000 -p 1e-5 -d 10600 -D 43100 -k 200 -gdir 0 1 0
 ```
##### For single-compartment domains
```bash
singularity exec -B $PWD writable_fenics_stable.simg python3 PreprocessingOneCompt.py -o onecompt_files.h5
singularity exec -B $PWD writable_fenics_stable.simg mpirun -n 6 python3 GCloudDmriSolver.py -f onecompt_files.h5 -M 0 -b 1000 -d 10600 -D 43100 -k 200 -gdir 1 0 0 
```
## Working with FEniCS-HPC Singularity Image
