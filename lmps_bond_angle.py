#!/usr/bin/env python

import sys

# inPOSCAR = 'POSCAR'
# inPOSCAR = sys.argv[1]

if len(sys.argv) < 2:
    inPOSCAR = 'POSCAR'
elif len(sys.argv) == 2:
    inPOSCAR = sys.argv[1]
else:
    print('ERROR: Too many command line arguments specified.')
    sys.exit()

out_data = inPOSCAR + '.bond_angle_info.txt'

def read_poscar(inPOSCAR):
    file_poscar = open(inPOSCAR, 'r')
    section = ''
    title = []
    multiplier = []
    cell_vecs = []
    coords = []
    water_coords = []
    
    line_count = 1
    atom_count = 1
    for line in file_poscar:
        if line_count == 1:
            section = 'title'
            title = line.rstrip('\n')
        if line_count == 2:
            section = 'multiplier'
            multiplier = line.strip('\n').split()
            x_vec = next(file_poscar).strip('\n').split()
            y_vec = next(file_poscar).strip('\n').split()
            z_vec = next(file_poscar).strip('\n').split()
            cell_vecs = [x_vec, y_vec, z_vec]
            line_count = 5
        if line_count == 6:
            elements = line.strip('\n').split()
        if line_count == 7:
            atoms = line.strip('\n').strip()
        if line_count == 8:
            if line.lower().startswith('s'):
                sel_dyn = 'Selective Dynamics'
            else:
                line_count += 1
        if line_count == 9:
            section = 'coord type'
            if line.lower().startswith('c'):
                coord_type = 'Cartesian'
            elif line.lower().startswith('d'):
                coord_type = 'Direct'
            else:
                sys.exit('ERROR: Incorrect coordinate type! Please change to Cartesian or Direct!')
        if line_count == 10:
            section = 'coords'
        line_count += 1


        if section == 'coords':
            coord = [atom_count] + line.strip('\n').split()
            if coord[-1] == '0':
                water_coords.append(coord)
            coords.append(coord)
            atom_count += 1
    
    file_poscar.close()
    
    return water_coords

def get_water_bonds_angles(water_coords):
    bonds_list = []
    angles_list = []
    
    if len(water_coords)%3 == 0:
        num_waters = len(water_coords)/3
    else:
        sys.exit('ERROR: The number of Oxygens and Hydrogens does not sum to a multiple of 3!')

    print(str(int(num_waters)) + ' water molecules were identified.')

    o_coords = water_coords[:int(num_waters)]
    h_coords = water_coords[int(num_waters):]
    if len(h_coords) != 2*len(o_coords):
        sys.exit('ERROR: There is an inconsistency between the number of Oxygens and Hydrogens!')
        
    for o_coord in o_coords:
        for h_coord in h_coords:
            if o_coord[4] == h_coord[4]:
                bonds_list.append([o_coord[0], h_coord[0]])
                
    print(str(len(bonds_list)) + ' bonds were identified.')
 
    for bond1 in bonds_list:
        for bond2 in bonds_list:
            if bond1[0] == bond2[0] and bond1[1] != bond2[1] and [bond2[1], bond1[0], bond1[1]] not in angles_list:
                angles_list.append([bond1[1], bond1[0], bond2[1]])
                
    print(str(len(angles_list)) + ' angles were identified.')

    with open(out_data, 'w') as f:
        f.write(' Bonds \n')
        f.write(' \n')
        bond_counter = 1
        for bond in bonds_list:
            f.write(str(bond_counter) + ' 1 ' + str(bond[0]) + ' ' + str(bond[1]) + '\n')
            bond_counter += 1
        f.write(' \n')
        f.write(' Angles \n')
        f.write(' \n')
        angle_counter = 1
        for angle in angles_list:
            f.write(str(angle_counter) + ' 1 ' + str(angle[0]) + ' ' + str(angle[1]) + ' ' + str(angle[2]) + '\n')
            angle_counter += 1
    
water_coords = read_poscar(inPOSCAR)
get_water_bonds_angles(water_coords)
