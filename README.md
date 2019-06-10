# Install Singularity
```bash
sudo apt-get singulariy-container
```
# Build images with Ubuntu
#### FEniCS
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_recipe_FEniCS_DMRI
sudo singularity build -w writable_fenics_dmri.simg Singularity_recipe_FEniCS_DMRI
```
#### FEniCS-HPC
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/build_fenics_hpc_image_recipe
sudo singularity build -w writable_fenics_hpc.simg build_fenics_hpc_image_recipe
```
### Test if mpi works correctly with the FEniCS-HPC image
```bash
wget https://raw.githubusercontent.com/wesleykendall/mpitutorial/gh-pages/tutorials/mpi-hello-world/code/mpi_hello_world.c
singularity exec -B $PWD writable_fenics_hpc_dmri.simg mpicc mpi_hello_world.c -o mpi_hello_world
mpirun -n 3 singularity exec -B $PWD writable_fenics_hpc.simg  mpi_hello_world
```
The results would be
```bash
Hello world from processor dmri, rank 0 out of 3 processors
Hello world from processor dmri, rank 1 out of 3 processors
Hello world from processor dmri, rank 2 out of 3 processors
```
