#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

log_filename = sys.argv[1]

equilibration_time = 2   # ns
log_write_freq = 1000    # write Lz to log file every 1000 timesteps

first_production_Lz = equilibration_time * 1000000 / log_write_freq + 1

def read_log(log_filename):
    log_file = open(log_filename,'r')
    section = ''
    subsection = ''
    
    nvt_file = open('nvt_data.txt','w+')
    nve_file = open('nve_data.txt','w+')
    npt_file = open('npt_data.txt','w+')
    
    for line in log_file:
        if line.startswith('------------ beginning equilibration (const Vol)------------'):
            section = 'nvt'
            continue
        elif line.startswith('------------ check energy conservation (only nve)------------'):
            section = 'nve'
            continue
        elif line.startswith('------------ beginning npt (use drag)------------'):
            section = 'npt'
            continue
        elif line.startswith('Loop'):
            section = ''
            subsection = ''
            continue

        if section == 'nvt':
            if line.startswith('Step'):
                subsection = 'nvt_data'
                line_data = line.split()
                nvt_temp_index = line_data.index('Temp')
                step = line_data[0]
                temp = line_data[nvt_temp_index]
                nvt_file.write('%-10s %-10s \n' % (step, temp))
                continue
        elif section == 'nve':
            if line.startswith('Step'):
                subsection = 'nve_data'
                line_data = line.split()
                nve_enrg_index = line_data.index('TotEng')
                step = line_data[0]
                enrg = line_data[nve_enrg_index]
                nve_file.write('%-10s %-10s \n' % (step, enrg))
                continue
        elif section == 'npt':
            if line.startswith('Step'):
                subsection = 'npt_data'
                line_data = line.split()
                npt_lz_index = line_data.index('Lz')
                step = line_data[0]
                lz = line_data[npt_lz_index]
                npt_file.write('%-10s %-10s \n' % (step, lz))
                continue

        if subsection == 'nvt_data':
            line_data = line.split()
            step = line_data[0]
            temp = line_data[nvt_temp_index]
            nvt_file.write('%-10s %-10s \n' % (step, temp))
        elif subsection == 'nve_data':
            line_data = line.split()
            step = line_data[0]
            enrg = line_data[nve_enrg_index]
            nve_file.write('%-10s %-10s \n' % (step, enrg))
        elif subsection == 'npt_data':
            line_data = line.split()
            step = line_data[0]
            lz   = line_data[npt_lz_index]
            npt_file.write('%-10s %-10s \n' % (step, lz))

    log_file.close()
    nvt_file.close()
    nve_file.close()
    npt_file.close()

read_log(log_filename)

plt.switch_backend('agg')

x = []
y = []
with open('nvt_data.txt') as f:
    header = f.readline()
    lines = f.readlines()
    x = [int(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
fig = plt.figure()
plt.plot(x, y)
plt.title('NVT Temp vs. Step')
plt.xlabel('Time Step')
plt.ylabel('Temperature [K]')
plt.savefig('nvt_plot.png')
plt.close(fig)

x = []
y = []
with open('nve_data.txt') as f:
    header = f.readline()
    lines = f.readlines()
    x = [int(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
fig = plt.figure()
plt.plot(x, y)
plt.title('NVE Energy vs. Step')
plt.xlabel('Time Step')
plt.ylabel('Energy [eV]')
plt.savefig('nve_plot.png')
plt.close(fig)

x = []
y = []
with open('npt_data.txt') as f:
    header = f.readline()
    lines = f.readlines()
    x = [int(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
fig = plt.figure()
plt.plot(x, y)
plt.title('NPT Lz vs. Step')
plt.xlabel('Time Step')
plt.ylabel('Length of z-axis [A]')
plt.savefig('npt_plot.png')
plt.close(fig)

with open('avg_lz.txt','w+') as f:
    f.write('The average Lz is: %.16f' % np.mean(np.array(y[first_production_Lz:])))
    f.write('The Lz st.dev. is: %.16f' % np.stdev(np.array(y[first_production_Lz:])))


