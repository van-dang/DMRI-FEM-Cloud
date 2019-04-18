# Copyright Van-Dang NGUYEN 2019
import os, sys

mshfile=sys.argv[1]

fext=mshfile[len(mshfile)-3:len(mshfile)]
if not(fext =="msh"):
  print("Sorry!!! The program supports .msh only! Please check the input file format...");
  exit(0);

print("Extracting cell markers from: "+mshfile+" ...");

filewithoutextention = os.path.splitext(os.path.basename(mshfile))[0];

f = open(mshfile, "r")
lineList = f.readlines()
f.close()

lastline = lineList[len(lineList)-2].split(" ")
last_id = int(lastline[0])

first_id = -100;

for x in lineList:
  x = x.split(" ")
  if len(x) == len(lastline):
        first_id = int(x[0])
        break;

if len(lastline)==8:
        dim = 2;
elif len(lastline)==9:
        dim = 3;
else:
        print("Invalid dimension. Please double check msh2xml!");

numcells = last_id - first_id + 1

outfile = 'volmk_'+filewithoutextention+'.xml'
filename = open(outfile,'w');
filename.write('<?xml version="1.0"?>\n')
filename.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
filename.write('  <mesh_function>\n')
filename.write('    <mesh_value_collection type="uint" dim="'+str(dim)+'" size="'+str(numcells)+'">\n');

for x in lineList:
  x = x.split(" ")
  if len(x) == len(lastline):
	  filename.write('      <value cell_index="'+str(int(x[0])-first_id)+'" local_entity="0" value="'+x[4]+'" />\n');

filename.write('    </mesh_value_collection>\n');
filename.write('  </mesh_function>\n');
filename.write('</dolfin>\n');
filename.close();

print("Extracted successfully on: "+str(numcells)+" elements")
print("Wrote to: "+outfile)