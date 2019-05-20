# Install Singularity
```bash
sudo apt-get singulariy-container
```
# Re-build images if neccessary
#### FEniCS
```bash
sudo singularity build fenics_stable.simg docker://fenicsproject/stable
```
#### FEniCS-HPC
```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/singularity_images/Singularity_DMRI_recipe
sudo singularity build fenics-hpc-dmri.simg Singularity_DMRI_recipe
```

# Download existing images
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics_stable.simg
```

# Change to writable mode
```bash
sudo singularity build --writable writable_fenics-hpc-dmri.simg fenics-hpc-dmri.simg
sudo singularity build --writable writable_fenics_stable.simg fenics_stable.simg
```

# Install neccessary packages to the existing image
```bash
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get install zip unzip gmsh
```
