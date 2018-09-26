#!/usr/bin/env python

#####
# 
# Cameron Bodenschatz
# Getman Research Group
# 26 September 2018
#
# USAGE: In LAMMPS job directory: lammps_frames.py output_file dump_file sort_method num_frames
#
# Acceptable sort_method values: energy, time
#
#####

import argparse, os, re
from operator import itemgetter

path = os.getcwd()

output_file = [f for f in os.listdir(path) if re.search(r'(.+).o([0-9]+)', f)][-1]
print('LAMMPS Output File: {}'.format(str(output_file)))

dump_file = [f for f in os.listdir(path) if re.search(r'dump.(.+).lammpstrj', f)][0]
print('LAMMPS Dump File: {}'.format(str(dump_file)))

parser = argparse.ArgumentParser()
parser.add_argument('--sort_method', '-s', default='time', choices=['time','energy'], help='Frame sort method.')
parser.add_argument('--num_frames', '-n', default=10, help='Number of frames.')
args = parser.parse_args()
sort_method = args.sort_method
num_frames = args.num_frames


def read_output(output_fn, production_start='2000000'):
    output_file = open(output_fn, 'r')
    section = ''
    subsection = ''

    production_data = []

    for line in output_file:

        if line.startswith('------------ beginning nvt ----------------------------------'):
            section = 'nvt'
            continue

        if line.startswith('Step'):
            subsection = 'thermo_data'
            continue

        if line.startswith(' '+production_start):
            subsection = 'production'
            continue

        if line.startswith('Loop'):
            subsection = 'job_details'
            continue

        if section is 'nvt' and subsection is 'production':
            if line.startswith('[node') or line.startswith(' [node'):
                continue
            else:
                row = line.rstrip('\n').split()
                production_data.append([row[0],row[3]])

    output_file.close()

    return production_data

prod_energies = read_output(output_file)       

def sort_energies(prod_energies, sort_method, num_frames):
    df = len(prod_energies)/num_frames
    if sort_method == 'time':
        frames = [];
        for i in range(0,num_frames):
            frames.append(prod_energies[int(df*(i+0.5)-1)][0])
    elif sort_method == 'energy':
        sort_energy = sorted(prod_energies, key=itemgetter(1))
        frames = [];
        for i in range(0,num_frames):
            frames.append(sort_energy[int(df*(i+0.5)-1)][0])
    else:
        print('ERROR: Incorrect sort method.')
        return
    return frames

frame_nums = sort_energies(prod_energies,sort_method,num_frames)

print('Creating files for the following trajectory timesteps: {}'.format(frame_nums))

# frame_nums = ['25000000','27500000','30000000','32500000','35000000','37500000','40000000','42500000','45000000','47500000']
#-----------------------------------------------------------------------------#

dump = [line.rstrip('\n') for line in open(dump_file)]

for i in range(0,len(dump)):
    for j in range(0,len(frame_nums)):
        if dump[i] == frame_nums[j]:
            
            f = open(str(j)+'_'+frame_nums[j], 'w+')
            f.writelines("%s\n" % dump[i-1])
            
            k = i;
            while dump[k] != 'ITEM: TIMESTEP':
                f.writelines("%s\n" % dump[k])
                k = k+1;
                if k >= len(dump):
                    break
            f.close()
            
del frame_nums
del dump

