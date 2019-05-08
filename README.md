# Buid images if neccessary
```bash
sudo singularity build fenics_stable.simg docker://fenicsproject/stable
sudo singularity build fenics-hpc-dmri.simg build_image_source_dmri
```

# Download images from the terminal
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics_stable.simg
```

# Change to writable mode
```bash
sudo singularity build --writable writable_fenics-hpc-dmri.simg fenics-hpc-dmri.simg
sudo singularity build --writable writable_fenics_stable.simg fenics_stable.simg
```

# Install packages to the image
```bash
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get update
sudo singularity exec --writable writable_fenics_stable.simg sudo apt-get install zip unzip gmsh
```
