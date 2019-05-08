# Compile
```bash
singularity exec -B $PWD ../../../writable_fenics-hpc-dmri.simg make -j 8
```
# Download the existing meshes
```bash
wget https://github.com/van-dang/NeuronVolumeMeshes/raw/master/pyramidals/05a_pyramidal8aACC.msh.zip
unzip 05a_pyramidal8aACC.msh.zip
```

# Convert .gmsh to .xml
```bash
wget https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.py
python python dolfin-convert.py 05a_pyramidal8aACC.msh 05a_pyramidal8aACC.xml
```

# Execute the demo
```bash
mpirun -n 8 singularity exec -B $PWD ../../../writable_fenics-hpc-dmri.simg ./demo -m 05a_pyramidal8aACC.xml 
```

```bash
wget https://raw.githubusercontent.com/van-dang/MRI-Cloud/fenics-hpc-solvers/build_image_source_dmri
sudo singularity build fenics-hpc-dmri.simg build_image_source_dmri
wget https://github.com/van-dang/MRI-Cloud/raw/fenics-hpc-solvers/fenics-hpc-dmri.simg
```
