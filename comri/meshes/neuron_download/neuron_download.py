spindle_list=[
'03a_spindle2aFI',
'03a_spindle6aFI',
'03b_spindle4aACC',
'03b_spindle5aACC',
'03b_spindle6aACC',
'03b_spindle7aACC',
'04b_spindle3aFI',
'05b_spindle5aFI',
'06b_spindle8aACC',
'07b_spindle9aACC',
'08a_spindle13aACC',
'09o_spindle7aFI',
'09o_spindle8aFI',
'10a_spindle18aACC',
'12a_spindle19aACC',
'12o_spindle9aFI',
'13o_spindle10aFI',
'15o_spindle12aFI',
'16o_spindle13aFI',
'19o_spindle14aFI',
'21o_spindle15aFI',
'23o_spindle16aFI',
'25o_spindle17aFI',
'26o_spindle18aFI',
'27o_spindle19aFI',
'28o_spindle20aFI',
'28o_spindle21aFI',
'29o_spindle22aFI',
'30o_spindle23aFI',
];

pyramidal_list=[
'02a_pyramidal2aFI',
'02b_pyramidal1aACC',
'02b_pyramidal1aFI',
'03a_pyramidal9aFI',
'03b_pyramidal2aACC',
'03b_pyramidal3aACC',
'03b_pyramidal3aFI',
'03b_pyramidal4aFI',
'03b_pyramidal9aFI',
'04a_pyramidal4aACC',
'04a_pyramidal5aACC',
'04b_pyramidal5aFI',
'04b_pyramidal6aACC',
'04b_pyramidal6aFI',
'04b_pyramidal7aACC',
'05a_pyramidal10aACC',
'05a_pyramidal8aACC',
'05b_pyramidal7aFI',
'05b_pyramidal8aFI',
'05b_pyramidal9aACC',
'06a_pyramidal11aACC',
'06b_pyramidal10aFI',
'06b_pyramidal12aACC',
'07a_pyramidal13aACC',
'07b_pyramidal14aACC',
'08o_pyramidal11aFI',
'10a_pyramidal15aACC',
'11a_pyramidal16aACC',
'11o_pyramidal12aFI',
'17o_pyramidal13aFI',
'18o_pyramidal14aFI',
'20o_pyramidal15aFI',
'22o_pyramidal16aFI',
'24o_pyramidal17aFI',
'25o_pyramidal18aFI',
'31o_pyramidal19aFI',
];

neuron_id = 1;
neuron_type = 'spindles';

if neuron_type == 'spindles':
    neuron_list = spindle_list;
if neuron_type == 'pyramidals':
    neuron_list = pyramidal_list;
    
neuron_name = neuron_list[neuron_id];
print(neuron_name)
neuron_dir='https://raw.githubusercontent.com/van-dang/NeuronVolumeMeshes/master/'+neuron_type+'/'+neuron_name+'.msh.zip'

import os

os.system("wget "+neuron_dir)
