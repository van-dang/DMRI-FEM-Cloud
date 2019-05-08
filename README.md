
```bash
sudo singularity build fenics_stable.simg docker://fenicsproject/stablesudo singularity build fenics-hpc-dmri.simg build_image_source_dmri
```

# Download images from the terminal
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics_stable.simg
```

# Change to writable mode
```bash
sudo singularity build --writable fenics-hpc-dmri-write.simg fenics-hpc-dmri.simg
sudo singularity build --writable fenics_stable-write.simg fenics_stable.simg
```
