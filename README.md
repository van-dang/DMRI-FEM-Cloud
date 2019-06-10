The following instructions were tested with Ubuntu 18.04
# Requirements:
* Singularity container
```bash
sudo apt-get singulariy-container
```
* OpenMPI which is compatible with the FEniCS-HPC image. Currently, it is openmpi-4.0.0

# Build the FEniCS-HPC image
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/build_fenics_hpc_image_recipe
sudo singularity build -w writable_fenics_hpc.simg build_fenics_hpc_image_recipe
```
# Test if mpi works correctly with the FEniCS-HPC image
```bash
wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
singularity exec -B $PWD writable_fenics_hpc.simg mpicc mpi_hello_world.c -o mpi_hello_world
mpirun -n 3 singularity exec -B $PWD writable_fenics_hpc.simg ./mpi_hello_world
```
The results would be
```bash
Hello world from processor dmri, rank 0 out of 3 processors
Hello world from processor dmri, rank 1 out of 3 processors
Hello world from processor dmri, rank 2 out of 3 processors
```
In case the openmpi is not compatible between the hosted machine and the image, it does not work correctly with a multi-node system. The following command can be used for one node.
```bash
singularity exec -B $PWD writable_fenics_hpc.simg mpirun -n 3 ./mpi_hello_world
```

