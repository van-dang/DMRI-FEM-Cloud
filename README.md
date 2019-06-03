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
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_recipe_FEniCS_HPC_DMRI
sudo singularity build -w writable_fenics_hpc_dmri.simg Singularity_recipe_FEniCS_HPC_DMRI
```
# Download existing images
#### FEniCS
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/writable_fenics_dmri.simg
```
#### FEniCS-HPC
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/writable_fenics_hpc_dmri.simg
```

Note: In addition to FEniCS/FEniCS-HPC platforms, these images were built with GMSH
