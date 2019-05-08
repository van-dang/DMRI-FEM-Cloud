# Install dependencies
```bash
sudo apt-get install singularity-container unzip mpich
```

# Download existing images
```bash
wget https://github.com/van-dang/MRI-Cloud/raw/singularity_images/fenics-hpc-dmri.simg
```

# Change to writable mode
```bash
sudo singularity build --writable writable_fenics-hpc-dmri.simg fenics-hpc-dmri.simg
```

# Install packages to the existing image
```bash
sudo singularity exec --writable writable_fenics-hpc-dmri.simg sudo apt-get update
sudo singularity exec --writable writable_fenics-hpc-dmri.simg sudo apt-get install zip unzip gmsh
```

# Download the solvers
```bash
wget https://github.com/van-dang/MRI-Cloud/archive/fenics-hpc-solvers.zip
```
# Compile
```bash
singularity exec -B $PWD ../../../writable_fenics-hpc-dmri.simg make -j 8
```
# Download the existing meshes
```bash
wget https://github.com/van-dang/NeuronVolumeMeshes/raw/master/pyramidals/04b_pyramidal7aACC.msh.zip
unzip 04b_pyramidal7aACC.msh.zip
```

# Convert .gmsh to .xml
```bash
wget https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.py
python dolfin-convert.py 04b_pyramidal7aACC.msh 04b_pyramidal7aACC.xml
```

# Execute the demo
```bash
singularity exec -B $PWD ../../../writable_fenics-hpc-dmri.simg mpirun -n 8 ./demo -m 04b_pyramidal7aACC.xml 
```

```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/fenics-hpc-solvers/build_image_source_dmri
sudo singularity build fenics-hpc-dmri.simg build_image_source_dmri
wget https://github.com/van-dang/MRI-Cloud/raw/fenics-hpc-solvers/fenics-hpc-dmri.simg
```
