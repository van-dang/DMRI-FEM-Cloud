# Install Singularity
```bash
sudo apt-get singulariy-container
```
# Build images with Ubuntu.
#### FEniCS
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_recipe_FEniCS_DMRI
sudo singularity build -w writable_fenics_stable.simg Singularity_recipe_FEniCS_DMRI
```
#### FEniCS-HPC
```bash
https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_recipe_FEniCS_HPC_DMRI
sudo singularity build -w writable_fenics-hpc-dmri.simg Singularity_recipe_FEniCS_HPC_DMRI
```

# Download existing images
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics_stable.simg
```

# Change to writable mode
```bash
sudo singularity build -w writable_fenics-hpc-dmri.simg fenics-hpc-dmri.simg
sudo singularity build -w writable_fenics_stable.simg fenics_stable.simg
```

# Install neccessary packages to the existing image
```bash
sudo singularity exec -w writable_fenics_stable.simg sudo apt-get install zip unzip gmsh
```
