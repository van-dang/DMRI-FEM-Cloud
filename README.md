# Simulation of diffusion MRI with Python notebooks

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

sudo apt-get install mpich singularity-container
```

## Create FEniCS Singularity Image

```bash
sudo singularity build --writable fenics_stable.simg docker://fenicsproject/stable
sudo singularity exec --writable fenics_stable.simg sudo apt-get update
sudo singularity exec --writable fenics_stable.simg sudo apt-get install zip unzip gmsh
```

## Copy Python solvers to the VM instance
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/PreprocessingOneCompt.py
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/PreprocessingMultiCompt.py
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/master/GCloudDmriSolver.py
```

## Execute the code with Singularity Image

```bash
singularity exec -B $PWD ../fenics_stable.simg python3 PreprocessingMultiCompt.py -o myfiles.h5
singularity exec -B $PWD ../fenics_stable.simg python3 GCloudDmriSolver.py -f myfiles.h5 -M 1 -b 1000 -k 200 -gdir 0 1 0
```
