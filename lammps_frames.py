#!/usr/bin/env python

#####
# 
# Cameron Bodenschatz
# Getman Research Group
# 26 September 2018
#
# USAGE: In LAMMPS job directory: python lammps_frames.py {-l log_file} {-d dump_file} {-n num_frames}
#
# NOTE: The larger your Lammps dump file is, the longer this script will take. It can take a couple of minutes.
#
#####

import argparse, os, re
from operator import itemgetter

path = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument('--num_frames', '-n', default=10, help='Number of frames.')
parser.add_argument('--log_file', '-l', default=None, help='Log file name.')
parser.add_argument('--dump_file', '-d', default=None, help='Dump file name.')

args = parser.parse_args()
num_frames = args.num_frames
log_file = args.log_file
dump_file = args.dump_file

if log_file is None:
    log_file = [f for f in os.listdir(path) if re.search(r'log.(.+)', f) and f != 'log.lammps'][-1]

if dump_file is None:
    dump_file = [f for f in os.listdir(path) if re.search(r'dump.(.+).lammpstrj', f)][0]

print('LAMMPS Output File: {}'.format(str(log_file)))
print('LAMMPS Dump File: {}'.format(str(dump_file)))


def read_output(log_fn, production_start='2000000'):
    log_file = open(log_fn, 'r')
    section = ''
    subsection = ''

    production_data = []

    for line in log_file:

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

    log_file.close()

    return production_data

prod_energies = read_output(log_file)       

def sort_energies(prod_energies, num_frames):
    df = len(prod_energies)/num_frames
    frames = [];
    for i in range(0,num_frames):
        frames.append(prod_energies[int(df*(i+0.5)-1)][0])
    return frames

frame_nums = sort_energies(prod_energies,num_frames)

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

