The following instructions were tested with Ubuntu 18.04
# Requirements:
* Singularity container
```bash
sudo apt-get singulariy-container
```

# Build the FEniCS-HPC image
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_recipe_FEniCS_HPC_DMRI
sudo singularity build -w writable_fenics_hpc.simg build_fenics_hpc_image_recipe
```
# Test if mpi works correctly with the FEniCS-HPC image
```bash
wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
singularity exec -B $PWD writable_fenics_hpc.simg mpicc mpi_hello_world.c -o mpi_hello_world
singularity exec -B $PWD writable_fenics_hpc.simg mpirun -n 3 ./mpi_hello_world
```
The results would be
```bash
Hello world from processor dmri, rank 0 out of 3 processors
Hello world from processor dmri, rank 1 out of 3 processors
Hello world from processor dmri, rank 2 out of 3 processors
```
# Note
For a multi-node system, openmpi needs to be compatible between the hosted machine and the image to launch with many processors beyond one node. It requires to install the same version of openmpi on the hosted machine and the image. The command to launch the demo is
```bash
mpirun -n 30 singularity exec -B $PWD writable_fenics_hpc.simg ./mpi_hello_world
```

