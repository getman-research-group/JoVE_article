#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

//Constants
#define eul 2.718281828
#define pi 3.14159265359
#define boltzmann 8.6173324E-5 // eV/K

#define nth_atom_limit 1 //Not equal to 1 only for random dihedral angles. For fixed dihedral angles, the limit is 1 attempt.
#define third_atom_limit 100 //Actual limit will be greater than 1
#define second_atom_limit 200 //Actual limit will be greater than 1 and greater than third_atom_limit
#define first_atom_limit 1000 //Actual limit will be greater than 1 and greater than second_atom_limit

#define initial_rand_seed -pi


/***********************************CHANGE DIRECTORY HERE************************************/

char file_prefix[1024] = "";

//LEAVE BLANK IF RUNNING ON PALMETTO
//The system automatically adds the pwd to any file names

/***********************************CHANGE DIRECTORY HERE************************************/








  ///////////////////////////////
 ////    Current Problems   ////
///////////////////////////////
/*
    1) The shared elements bug is still there. It has to be somewhere in the replacement_sorter() function, but I haven't figured it out yet.

    2) The translate and rotate functionalities don't work at all with large (>3 atom) molecules.

    3) The bond angle perturbation routine isn't finished yet.
*/





  /////////////////////////////////////
 ////    Structure definitions    ////
/////////////////////////////////////

//POSCAR

typedef struct
{
    double scale_factor;
    double latt_vec[9];
    int num_types; //number of different types of atoms
    int* par_num; //contains the number of each type of atom
                  //space will be allocated on an individual basis
    int par_count;
    char par_ident[112][3]; //contains atom symbols
    double* par_pos; //contains coordinates of all particles in all instances of the molecule.
                     //space will be allocated on an individual basis (I forgot what this means. I'm not sure it's true, anyway.)
                        //for each instance of the molecule. Atoms WILL NOT be sorted according to identity.
    double energy;

    int rm_terminated; /**Replacement**/

    int move_type; // 0 = translation, 1 = rotation, 2 = bond angle perturbation

                     //atoms will be ordered as they are in the .fh file, and this order will be repeated
    int* molecular_ID; //contains the molecular ID tag of every atom, so that molecular bonding information is preserved. # members = num_atoms, and the atom order in this array is the same as the atom order in the position array

    int* molecular_type_tag; //contains the molecular type tag for every atom. This allows the program to determine which atoms belong to which molecular type (e.g., CO2 vs H20) upon reading a resubmitted POSCAR file

    int* atom_ID;
    int* bonded_to;
    int* angle;

} POSCAR_struct;


//Placed particles

typedef struct
{
    int num_atoms;
    char atom_ident[112][3]; //contains atomic symbol for each atom in the molecule in order of placement

    //Relative location information//

    int* bonded_to; //contains the ID number (placement order, e.g. 1 (first), 2 (second), etc.) of each atom
    double* lengths_bond; //contains the bond length, in Angstroms, between every adjacent pair of atoms.
    int* third_atom; //contains the ID number of each atom used as the third point in a bond angle
    double* angles_bond; //contains the bond angles for each adjacent trio of atoms
    int* fourth_atom; //contains the ID number of each atom used as the fourth point for a dihedral angle
    double* dihedral; //contains the dihedral angle between every set of four adjacent atoms

    double* bond_angle_ranges;

} molecule_struct;


typedef struct
{
    int max_number; //maximum number of placed molecules of type X
    int actual_number_placed; //actual_number_placed may be less than max_number because of space limitations

    int molecular_type_tag; //This keeps track of the type of molecule in the POSCAR (and possibly later XDATCAR) files.
                            //It'll be used to simplify the RESUBMIT routine

    molecule_struct* structure; //pointer to the structure containing the relative location information of a particular molecule

    double* min_distance_array;

    int num_atom_types; //number of different types of atoms in the molecule
    int* num_each_type; //number of each type of atom in the molecule
        //combine these two with max_num to calculate num_atoms
    int* max_num_each_type;
    char type_ident[112][3]; //contains atomic symbol for each atom type in the molecule in order of placement

    int close_enough; //==0 if current atom not close to surface atom or another placed atom (outside of its own molecule)
                        // ==1 if the current atom is close enough to another atom outside of its own molecule

    int* atom_ID;
    int* bonded_to; //stores bonding information for LAMMPS file
    int* angle; //stores angle information for LAMMPS file

    double* coords; //contains the coordinates of all the atoms of the first molecule, followed by all the atoms of the second instance of that molecule, followed by the third instance...
                    //The coordinates are ordered in groups of three by the order in which they are placed in the placement routines.
                        //e.g., if CH30H is placed, and the atomic placement order is C,O,H,H,H,H, that is the order in which the coordinates will appear.
                    //The above is only true if a RESUBMIT_POSCAR is NOT being used. This is only true if initial particle placement takes place
                    //Otherwise, the atoms are grouped by molecule and subgrouped by atom type
///////////////////////////////////////////////////////////////////////////////
    int rotate_about_this_atom; /**Replacement**/ //Uses index numbers (starts at 0) to identify which atom in the replaced molecule will be the point about which the rest of the molecule is rotated

    int translate; /**Replacement**/ // These are set equal to 1 or 0 to indicate the type of random move being made
    int rotate; /**Replacement**/ // ^^
    int BAP; /**Replacement**/ // ^^

    int* molecular_ID; //contains the molecular ID tag of every atom in this particular molecule type, so that molecular bonding information is preserved
                            //I'll need to increment the molecular_ID of atoms in molecules beyond the first one by the number of molecules of all previous types already placed (molecule_count in the placed_particles structure)

    int* replaced_par_pos_indices; /**Replacement**/ //Contains the position array indices of all atoms in the molecule being replaced
    int replacing; /**Replacement**/ //set to 1 to indicate molecular_ID[] shouldn't be filled starting at 0.
    int replaced_molecular_ID; /**Replacement**/ //the molecular ID of the molecule being replaced

} molecule_coords;


typedef struct
{
    int num_types_molecules; //the number of different kinds of molecules in the simulation

    int master_num_types_atoms;
    int master_max_num_each_type[112];
    int* master_actual_num_each_type;
    char master_ident[112][3]; //contains POSCAR atom identities, too, but these are placed here later

    int molecule_count;

    double rand_seed; //needed to ensure PRNG doesn't start with the same seed. This stores the last random value generated and uses it as the next seed.

    molecule_coords* molecule_info; //contains pointers to the molecule_coords structure for each different type of molecule
                                    //must allocate space for the number of POINTERS needed

    int terminate_replacing; /**Replacement**/ //set to 1 to end replacement routine once a placement routine failure occurs

} placed_particles;


int num_elements = 96;
//H: 0.31, 1.10
double covalent_radii[96] ={0.18,  //0.31
                            0.28,
                            1.28,
                            0.96,
                            0.84,
                            0.76,
                            0.71,
                            1.42,  //0.66
                            0.57,
                            0.58,

                            1.66,
                            1.41,
                            1.21,
                            1.11,
                            1.07,
                            1.05,
                            1.02,
                            1.06,
                            2.03,
                            1.76,

                            1.70,
                            1.60,
                            1.53,
                            1.39,
                            1.39,
                            1.32,
                            1.22,
                            1.22,
                            1.20,
                            1.19,

                            1.20,
                            1.20,
                            1.16,
                            2.20,
                            1.95,
                            1.80,
                            1.75,
                            1.64,
                            1.54,
                            1.47,

                            1.46,
                            1.42,
                            1.39,
                            1.45,
                            1.44,
                            1.42,
                            1.39,
                            1.39,
                            1.38,
                            1.39,

                            1.40,
                            2.44,
                            2.15,
                            1.87,
                            1.75,
                            1.70,
                            1.62,
                            1.51,
                            1.44,
                            1.41,

                            1.36,
                            1.36,
                            1.32,
                            1.45,
                            1.46,
                            1.48,
                            1.40,
                            1.50,
                            1.50,
                            2.60,

                            2.21,
                            2.07,
                            2.04,
                            2.03,
                            2.01,
                            1.99,
                            1.98,
                            1.98,
                            1.96,
                            1.94,

                            1.92,
                            1.92,
                            1.89,
                            1.90,
                            1.87,
                            2.15,
                            2.06,
                            2.00,
                            1.96,
                            1.90,

                            2.00,
                            1.96,
                            1.90,
                            1.87,
                            1.80,
                            1.69};

double vdw_radii[96] = {1.10,
                        1.40,
                        1.81,
                        1.53,
                        1.92,
                        1.70,
                        1.55,
                        1.52,
                        1.47,
                        1.54,

                        2.27,
                        1.73,
                        1.84,
                        2.10,
                        1.80,
                        1.80,
                        1.75,
                        1.88,
                        2.75,
                        2.31,

                        2.16,
                        1.87,
                        1.79,
                        1.89,
                        1.97,
                        1.94,
                        1.92,
                        1.84,
                        1.86,
                        2.10,

                        1.87,
                        2.11,
                        1.85,
                        1.90,
                        1.83,
                        2.02,
                        3.03,
                        2.49,
                        2.19,
                        1.86,

                        2.07,
                        2.09,
                        2.09,
                        2.07,
                        1.95,
                        2.02,
                        2.03,
                        2.30,
                        1.93,
                        2.17,

                        2.06,
                        2.06,
                        1.98,
                        2.16,
                        3.43,
                        2.68,
                        2.21,
                        2.12,
                        2.17,
                        2.10,

                        2.17,
                        2.16,
                        2.02,
                        2.09,
                        2.17,
                        2.09,
                        1.96,
                        2.02,
                        2.07,
                        1.97,

                        2.02,
                        2.20,
                        3.48,
                        2.83,
                        2.40,
                        2.35,
                        2.39,
                        2.29,
                        2.36,
                        2.29,

                        2.33,
                        2.37,
                        2.21,
                        2.29,
                        2.16,
                        2.35,
                        2.27,
                        2.42,
                        2.60,
                        2.37,

                        2.43,
                        2.40,
                        2.21,
                        2.43,
                        2.44,
                        2.45};

char* element_ID[96] ={"H",
                     "He",
                     "Li",
                     "Be",
                     "B",
                     "C",
                     "N",
                     "O",
                     "F",
                     "Ne",

                     "Na",
                     "Mg",
                     "Al",
                     "Si",
                     "P",
                     "S",
                     "Cl",
                     "Ar",
                     "K",
                     "Ca",

                     "Sc",
                     "Ti",
                     "V",
                     "Cr",
                     "Mn",
                     "Fe",
                     "Co",
                     "Ni",
                     "Cu",
                     "Zn",

                     "Ga",
                     "Ge",
                     "As",
                     "Se",
                     "Br",
                     "Kr",
                     "Rb",
                     "Sr",
                     "Y",
                     "Zr",

                     "Nb",
                     "Mo",
                     "Tc",
                     "Ru",
                     "Rh",
                     "Pd",
                     "Ag",
                     "Cd",
                     "In",
                     "Sn",

                     "Sb",
                     "Te",
                     "I",
                     "Xe",
                     "Cs",
                     "Ba",
                     "Lu",
                     "Hf",
                     "Ta",
                     "W",

                     "Re",
                     "Os",
                     "Ir",
                     "Pt",
                     "Au",
                     "Hg",
                     "Tl",
                     "Pb",
                     "Bi",
                     "Po",

                     "At",
                     "Rn",
                     "Fr",
                     "Ra",
                     "La",
                     "Ce",
                     "Pr",
                     "Nd",
                     "Pm",
                     "Sm",

                     "Eu",
                     "Gd",
                     "Tb",
                     "Dy",
                     "Ho",
                     "Er",
                     "Tm",
                     "Yb",
                     "Ac",
                     "Th",

                     "Pa",
                     "U",
                     "Np",
                     "Pu",
                     "Am",
                     "Cm"};




/*

Array_pos+1) element_ID covalent_radii[Angstroms] -- vdw_radii[Angstroms]

1) H   0.31 -- 1.10
2) He  0.28 -- 1.40
3) Li  1.28 -- 1.81
4) Be  0.96 -- 1.53
5) B   0.84 -- 1.92
6) C   0.76 -- 1.70
7) N   0.71 -- 1.55
8) O   0.66 -- 1.52
9) F   0.57 -- 1.47
10) Ne 0.58 -- 1.54

11) Na 1.66 -- 2.27
12) Mg 1.41 -- 1.73
13) Al 1.21 -- 1.84
14) Si 1.11 -- 2.10
15) P  1.07 -- 1.80
16) S  1.05 -- 1.80
17) Cl 1.02 -- 1.75
18) Ar 1.06 -- 1.88
19) K  2.03 -- 2.75
20) Ca 1.76 -- 2.31

21) Sc 1.70 -- 2.16
22) Ti 1.60 -- 1.87
23) V  1.53 -- 1.79
24) Cr 1.39 -- 1.89
25) Mn 1.39 -- 1.97
26) Fe 1.32 -- 1.94
27) Co 1.22 -- 1.92
28) Ni 1.22 -- 1.84
29) Cu 1.20 -- 1.86
30) Zn 1.19 -- 2.10

31) Ga 1.20 -- 1.87
32) Ge 1.20 -- 2.11
33) As 1.16 -- 1.85
34) Se 2.20 -- 1.90
35) Br 1.95 -- 1.83
36) Kr 1.80 -- 2.02
37) Rb 1.75 -- 3.03
38) Sr 1.64 -- 2.49
39) Y  1.54 -- 2.19
40) Zr 1.47 -- 1.86

41) Nb 1.46 -- 2.07
42) Mo 1.42 -- 2.09
43) Tc 1.39 -- 2.09
44) Ru 1.45 -- 2.07
45) Rh 1.44 -- 1.95
46) Pd 1.42 -- 2.02
47) Ag 1.39 -- 2.03
48) Cd 1.39 -- 2.30
49) In 1.38 -- 1.93
50) Sn 1.39 -- 2.17

51) Sb 1.40 -- 2.06
52) Te 2.44 -- 2.06
53) I  2.15 -- 1.98
54) Xe 1.87 -- 2.16
55) Cs 1.75 -- 3.43
56) Ba 1.70 -- 2.68
57) Lu 1.62 -- 2.21
58) Hf 1.51 -- 2.12
59) Ta 1.44 -- 2.17
60) W  1.41 -- 2.10

61) Re 1.36 -- 2.17
62) Os 1.36 -- 2.16
63) Ir 1.32 -- 2.02
64) Pt 1.45 -- 2.09
65) Au 1.46 -- 2.17
66) Hg 1.48 -- 2.09
67) Tl 1.40 -- 1.96
68) Pb 1.50 -- 2.02
69) Bi 1.50 -- 2.07
70) Po 2.60 -- 1.97

71) At 2.21 -- 2.02
72) Rn 2.07 -- 2.20
73) Fr 2.04 -- 3.48
74) Ra 2.03 -- 2.83
75) La 2.01 -- 2.40
76) Ce 1.99 -- 2.35
77) Pr 1.98 -- 2.39
78) Nd 1.98 -- 2.29
79) Pm 1.96 -- 2.36
80) Sm 1.94 -- 2.29

81) Eu 1.92 -- 2.33
82) Gd 1.92 -- 2.37
83) Tb 1.89 -- 2.21
84) Dy 1.90 -- 2.29
85) Ho 1.87 -- 2.16
86) Er 2.15 -- 2.35
87) Tm 2.06 -- 2.27
88) Yb 2.00 -- 2.42
89) Ac 1.96 -- 2.60
90) Th 1.90 -- 2.37

91) Pa 2.00 -- 2.43
92) U  1.96 -- 2.40
93) Np 1.90 -- 2.21
94) Pu 1.87 -- 2.43
95) Am 1.80 -- 2.44
96) Cm 1.69 -- 2.45


Sources:
1) http://webelements.com/periodicity/van_der_waals_radius/
2) http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_%28data_page%29
3) http://periodic.lanl.gov/list.shtml

**The d-block elements' values may be uncertain**

*/


  ////////////////////////////////////
 ////    Function Definitions    ////
////////////////////////////////////


/**THE FUNCTIONS APPEAR IN THE BODY OF THE PROGRAM IN THE FOLLOWING ORDER,
    AND ARE GROUPED BY THE HEADERS USED IN THIS SECTION**/

//Read input files
    int read_poscar(POSCAR_struct* input,
                    char* path,
                    int resubmit);

    placed_particles* read_internal_coords(FILE* master_input, char* file_prefix);



//Particle placement functions

    double* box_max(double sp_pos_carte[],
                    int sp_count,
                    POSCAR_struct* input);

    double* box_min(double sp_pos_carte[],
                    int sp_count,
                    POSCAR_struct* input);

    int particle_placement_shell(placed_particles* particle_list,
                                 POSCAR_struct* input,
                                 int zone,
                                 double min_scale,
                                 double max_scale,
                                 double z_min,
                                 double z_max,
                                 double maximum_translation);

    double* atom_placer(POSCAR_struct* input,
                        double* coords_pointer,
                        placed_particles* par_list,
                        molecule_coords* molec_info,
                        int atom_number,
                        int bonded_to,
                        double lengths_bond,
                        int third_atom,
                        double angles_bond,
                        int fourth_atom,
                        double dihedral,
                        double max_dim[],
                        double min_dim[],
                        double min_scale,
                        double max_scale,
                        int z_min,
                        int z_max);

    double* place_atom(placed_particles* par_list,
                       double max_dim[],
                       double min_dim[]);

    double* place_with_bond_length(placed_particles* par_list,
                                   double coords[],
                                   int bonded_to,
                                   double bond_length);

    double* place_with_bond_angle(placed_particles* par_list,
                                  double coords[],
                                  int bonded_to,
                                  double bond_length,
                                  int third_atom,
                                  double bond_angle);

    double* place_with_dihedral(placed_particles* par_list,
                                double coords[],
                                int bonded_to,
                                double bond_length,
                                int third_atom,
                                double bond_angle,
                                int fourth_atom,
                                double dihedral);

    double placement_checker_POSCAR(POSCAR_struct* POSCAR,
                                    molecule_coords* molec_info,
                                    int atom_number,
                                    double min_scale,
                                    double max_scale,
                                    double position_and_info[],
                                    int z_min,
                                    int z_max);


    double placement_checker_placedparticles(POSCAR_struct* POSCAR,
                                             placed_particles* par_list,
                                             molecule_coords* molec_info,
                                             int atom_number,
                                             double min_scale,
                                             double max_scale,
                                             double position_and_info[]);

    void fill_actual_number_placed(placed_particles* par_list, POSCAR_struct* resubmission);



//Output and organizing functions
    int organize_atom_idents(placed_particles* particle_list,
                             POSCAR_struct* input);

    void POSCAR_element_sorter(POSCAR_struct* POSCAR,
                               placed_particles* par_list);


    void output_sorter(POSCAR_struct* input,
                       placed_particles* par_list,
                       POSCAR_struct* output);


    void write_POSCAR(POSCAR_struct* output,
                      char* file_path,
                      placed_particles* par_list,
                      int finished);


    int return_radius_index(char idents[3]);

    void write_XDATCAR(char* file_prefix,
                       int num_POSCARS,
                       placed_particles* master,
                       int iteration);

    void write_LAMMPS(POSCAR_struct* input,
                      placed_particles* par_list,
                      int num_original_POSCAR_atom_types,
                      char file_extension[]);


//MC looping structure functions

    void energy_calc_simple(POSCAR_struct* surface,
                            int iteration_num,
                            int falling_particles);

    POSCAR_struct* random_move(POSCAR_struct* first_output,
                     POSCAR_struct* second_output,
                     placed_particles* replacement,
                     placed_particles* particle_list,
                     double prob_of_BAP,
                     double prob_of_rota,
                     double prob_of_trans,
                     double min_scale,
                     double max_scale,
                     double z_min,
                     double z_max,
                     int iteration_num,
                     int zone,
                     double maximum_translation,
                     int rm_termination_limit);

    void replacement_sorter(placed_particles* replacement,
                            POSCAR_struct* second_output,
                            POSCAR_struct* written,
                            int last_molecular_ID,
                            int selection);



//Helper functions
    void copyString(char* to,
                    char* from);

    double dist(double x1, double y1, double z1,
                double x2, double y2, double z2);

    double* num_gen(int runs,
                    placed_particles* par_list);

    double dot(double* vector1,
               double* vector2);

    double* create_unit_vector(double* head_pos,
                               double* tail_pos);

    int isEqual(int x,
                double y);

    void copy_POSCAR(POSCAR_struct* dest,
                     POSCAR_struct* source);

    double* cross(double a_pos[3],
                  double b_pos[3]);

    void unitize(double vector[3]);

    void VASP_energies(POSCAR_struct* blah,
                       int iteration_num,
                       int move_type,
                       int accepted,
                       int accepted_move_number,
                       int attempted_move_number,
                       int total_num_translate,
                       int total_num_rotate,
                       int finished,
                       double probability_of_rotation);

    int find_max_int(int array[], int num_members);




//Coordinate conversion functions
    double* frac_to_cartesian(double latt_vec[],
                              double pos_frac[],
                              int p_count,
                              double scale_factor);

    double* cartesian_to_frac(double latt_vec[],
                              double pos_cart[],
                              int p_count,
                              double scale_factor);


//Memory freeing functions

    void free_placed_particles(placed_particles* par_list);

    void free_POSCAR_struct(POSCAR_struct* poscar_file);


  //////////////////////////////
 ////    Main Code Block   ////
//////////////////////////////


int main()
{
   /*Variable, Pointer, and File Initialization*/
    int falling_particles = 1;
    int length = 0;
    int zone = 3;
    int num_atoms;

    int a, b; //counter variables

    char POSCAR_name[50]; // = (char*) malloc(50*sizeof(char));
    char* file_path = NULL;

    placed_particles* particle_list = NULL;

    POSCAR_struct* first_output = (POSCAR_struct*) malloc(sizeof(POSCAR_struct));

    FILE* master_input = NULL;

    char* master_input_file_path = NULL;

        length = strlen("master_input.txt") + strlen(file_prefix)+1;

        master_input_file_path = (char*) malloc(length*sizeof(char));

        sprintf(master_input_file_path, "%s%s", file_prefix, "master_input.txt");

        printf("\nmaster_input_file_path = %s", master_input_file_path);

    master_input = fopen(master_input_file_path,"r");

    double* dub = NULL;

    char* ident_holder;
    int* bonded_to;
    double* lengths_bond;
    double* third_atom;
    double* angles_bond;
    double* fourth_atom;
    double* dihedral;

    double z_min = 0;
    double z_max = 0;

    double max_scale = 1;
    double min_scale = 1;

    double trans_rotat_ratio = 0.5; //default value, changed soon.
    double prob_of_trans = 0.33; //default value
    double prob_of_rota = 0.33; //default value
    double prob_of_BAP = 0.33; //default value


    double temperature; //Kelvin
    int max_iterations = 100;
    double maximum_translation;


    //Variables for energies.txt file

    int move_type; // 0 = translation, 1 = rotation
    int accepted; // 1 = accepted, 0 = rejected
    int accepted_move_number = 1; // counter variable to keep track of the current number of moves that have been accepted by the MC routine
    int attempted_move_number = 0; // counter variable to keep track of the total number of moves attempted, both accepted and rejected (should be able to use the variable 'iteration')
    int total_num_translate = 0; // keeps track of the total number of SUCCESSFULL translation moves
    int total_num_rotate = 0; // keeps track of the total number of SUCCESSFULL rotation moves

    int resubmit;

    int rm_termination_limit = 50;
    int anneal = 0;
    double Tmax, Tmin;
    int anneal_begin = 0;
    int anneal_end = max_iterations;

    int random_waters = 0;

    char LAMMPS[64];


///////////////////////////////////////////////////////////////

    //Fills Initial POSCAR File

    if(master_input==NULL)
    {
        printf("\n\nCRASH: Unable to open master_input.txt\nMake sure the file path and name are correct\n\n");
        exit(1);
    }

    fscanf(master_input,"Falling particles (1 or 0): %i", &falling_particles);

    if((falling_particles!=0) && (falling_particles!=1))
    {
        printf("\nCRASH: Invalid value for falling_particles variable. Check value in master_input.txt.\n");
        exit(915);
    }

    fscanf(master_input,"\nResubmit (1 or 0): %i", &resubmit);

    if((resubmit!=0) && (resubmit!=1))
    {
        printf("\nCRASH: Invalid value for resubmit variable. Check value in master_input.txt.\n");
        exit(923);
    }

    fscanf(master_input,"\nWrite LAMMPS file (1 or 0): %i", &random_waters);

    if((random_waters!=0) && (random_waters!=1))
    {
        printf("\nCRASH: Invalid value for LAMMPS file variable. Check value in master_input.txt.\n");
        exit(938);
    }


    fscanf(master_input,"\nLAMMPS file extension: %s", LAMMPS);

    fscanf(master_input,"\nInput POSCAR name: %s", POSCAR_name);

    POSCAR_struct* Initial = (POSCAR_struct*) malloc(sizeof(POSCAR_struct));
    //RESUBMIT_POSCAR* Initial_resubmit = (RESUBMIT_POSCAR*) malloc(sizeof(RESUBMIT_POSCAR));


    fscanf(master_input,"\nNumber of different types of surface atoms: %i", &Initial->num_types);



    //printf("\nInitial->num_types = %i\n",Initial->num_types);

    length+=(strlen(POSCAR_name) + strlen(file_prefix)+1);

    file_path = (char*) malloc(length*sizeof(char));

    sprintf(file_path, "%s%s", file_prefix, &POSCAR_name);

//////////////////////////////////////////////////////////////////////

   /*****Call to read_poscar()*****/

    read_poscar(Initial, file_path, resubmit);  //fills Initial with initial input POSCAR information


    fscanf(master_input, "\nMove attempt termination limit: %i", &rm_termination_limit);

    if(rm_termination_limit<0)
    {
        printf("\nCRASH: Invalid value for move attempt termination limit. Check value in master_input.txt.\n");
        exit(955);
    }



    fscanf(master_input, "\nMaximum scaling factor: %lf", &max_scale);

    if(max_scale<0)
    {
        printf("\nCRASH: Invalid value for maximum scaling factor. Check value in master_input.txt.\n");
        exit(963);
    }


    fscanf(master_input, "\nMinimum scaling factor: %lf", &min_scale);

    if(min_scale<0)
    {
        printf("\nCRASH: Invalid value for minimum scaling factor. Check value in master_input.txt.\n");
    }

    fscanf(master_input, "\nMinimum z-coordinate: %lf", &z_min);

    fscanf(master_input, "\nMaximum z-coordinate: %lf", &z_max);

    fscanf(master_input, "\nMove probabilities [%%]\nTranslation: %lf", &prob_of_trans);

    fscanf(master_input, "\nRotation: %lf", &prob_of_rota);

    fscanf(master_input, "\nBond angle perturbation: %lf", &prob_of_BAP);

    if(((prob_of_BAP+prob_of_rota+prob_of_trans)>100) || ((prob_of_BAP+prob_of_rota+prob_of_trans)<0))
    {
        printf("\nInvalid values for move probabilities\nMake sure they sum to be less than 100 and greater than 0\nExiting");
        exit(1026);
    }



    fscanf(master_input, "\nTemperature [K]: %lf", &temperature);
    Tmin = temperature;

    if(temperature<0)
    {
        printf("\nCRASH: Invalid value for temperature. Check value in master_input.txt.\n");
        exit(991);
    }

    fscanf(master_input, "\nAnneal (1 or 0): %i", &anneal);

    if((anneal!=0) && (anneal!=1))
    {
        printf("\nCRASH: Invalid value for anneal variable. Check value in master_input.txt.\n");
        exit(991);
    }

    fscanf(master_input, "\nInitial temperature [K]: %lf", &Tmax);

    if(anneal==1)
    {
        temperature = Tmax;
    }

    if((Tmax<0) && (anneal==1))
    {
        printf("\nCRASH: Invalid value for initial annealing temperature. Check value in master_input.txt.\n");
        exit(1012);
    }

    fscanf(master_input, "\nBegin at iteration: %i", &anneal_begin);

    if((anneal_begin<0) && (anneal==1))
    {
        printf("\nCRASH: Invalid value for first annealing iteration. Check value in master_input.txt.\n");
        exit(1020);
    }

    fscanf(master_input, "\nEnd at iteration: %i", &anneal_end);

    if(((anneal_end<0) || (anneal_end<anneal_begin)) && (anneal==1))
    {
        printf("\nCRASH: Invalid value for last annealing iteration. Check value in master_input.txt.\n");
        exit(1028);
    }




    fscanf(master_input, "\nMax iterations: %i", &max_iterations);

    if(max_iterations<0)
    {
        printf("\nCRASH: Invalid value for maximum number of iterations. Check value in master_input.txt.\n");
        exit(1036);
    }


    fscanf(master_input,"\nMaximum translation distance [A]: %lf", &maximum_translation);

    if(maximum_translation<0)
    {
        printf("\nCRASH: Invalid value for maximum translation distance. Check value in master_input.txt.\n");
        exit(963);
    }


//Reads in information about placed particles

    particle_list = read_internal_coords(master_input, file_prefix);  //reads all necessary pdb file data into program

    particle_list->rand_seed = initial_rand_seed;


////////////////////////////////////////////////////////////////////////////////////////

//Fills in the remainder of the placed_particles struct in preparation for particle placement

    organize_atom_idents(particle_list, Initial);

//Begin particle placement routine
printf("\nBeginning fluid molecule placement routine");
    if(resubmit==0)
    {
        particle_placement_shell(particle_list, Initial, zone, min_scale, max_scale, z_min, z_max, maximum_translation);
    }
    else if(resubmit==1)
    {
        fill_actual_number_placed(particle_list, Initial);
    }


//Sorting program data and writing to a file (first POSCAR file)

    length = strlen(file_prefix) + strlen("POSCAR_0.POSCAR") + 1;

    char* first_POSCAR_file_path = (char*) malloc(length*sizeof(char));

    sprintf(first_POSCAR_file_path, "%s%s", file_prefix, "POSCAR_0.POSCAR");

    if(resubmit==0)
    {
        POSCAR_element_sorter(Initial, particle_list);
    }

    if(resubmit==0)
    {
        output_sorter(Initial, particle_list, first_output);
    }

    if(resubmit==1)
    {
        copy_POSCAR(first_output, Initial);
    }

    write_POSCAR(first_output, first_POSCAR_file_path, particle_list, 0);



    if(random_waters==1)
    {
        write_LAMMPS(first_output, particle_list, Initial->num_types, LAMMPS);
    }


///Beginning MC looping structure///
    int aa, bb, cc;

//Making new placed_particles structure for particle replacement
    placed_particles* replacement = (placed_particles*) malloc(sizeof(placed_particles));
    replacement->terminate_replacing = 0;

    energy_calc_simple(first_output, 0, falling_particles);

    VASP_energies(first_output, 0, 0, 1, 1, 1, 0, 0, 0, 0);


//Creating and allocating memory for the second_output POSCAR
    POSCAR_struct* second_output = (POSCAR_struct*) malloc(sizeof(POSCAR_struct));

    copy_POSCAR(second_output, first_output);


    POSCAR_struct* correct_output = NULL;

    length = strlen(file_prefix) + strlen("POSCAR_") + strlen(".POSCAR") + 1;

    char* POSCAR_file_path = (char*) malloc((length+sizeof(int)+sizeof(char)+1)*sizeof(char)); // (char*) malloc((strlen("D://POSCAR") + strlen(".POSCAR") + sizeof(int))*sizeof(char));

    int rm_termination_counter = 0;

    if(POSCAR_file_path==NULL)
    {
        printf("\nFailed to allocate POSCAR_file_path @ 1069!");
        exit(959);
    }

    char energies_path[1000];

    sprintf(energies_path, "%senergies.txt", file_prefix);

    //if(falling_particles==0)
    //{
    FILE* energies = NULL; //These three lines have the effect of clearing previous runs' data in preparation for the current run.

    energies = fopen(energies_path, "w");

    fclose(energies);
    //}

//Beginning main loop of MC
    int a_finished_on = 0;

    int attempt = 1;

    double* num_for_energy_decision;

    double ratio;

    int success = 0;

    write_XDATCAR(file_prefix, a_finished_on, particle_list, 0); //begins writing the XDATCAR file

    //printf("\nA crash after random water placement is OK.\n");

    for(a=0; a<max_iterations; a++)
    { //printf("\nreplacement->terminate_replacing = %i", replacement->terminate_replacing); printf("\tparticle_list->molecule_count = %i", particle_list->molecule_count);
        if((replacement->terminate_replacing==0) && (particle_list->molecule_count>0))
        { printf("\n\n\n\t\tATTEMPTING POSCAR_%i", a+1);
            correct_output = random_move(first_output, second_output, replacement, particle_list, prob_of_BAP, prob_of_rota, prob_of_trans, min_scale, max_scale, z_min, z_max, a, zone, maximum_translation, rm_termination_limit); //randomly moves one of the fluid molecules

            sprintf(POSCAR_file_path, "%s%s%i%s", file_prefix, "POSCAR_", a+1, ".POSCAR"); //creates POSCAR file name for future use
            write_POSCAR(correct_output, POSCAR_file_path, replacement, 0); //this state may or may not be accepted, but will eventually be overwritten by a state that is guaranteed to be accepted
                                                            //writes the current POSCAR
            particle_list->rand_seed = replacement->rand_seed; //transfers the current random seed from one structure to the current structure so that the PRNG doesn't begin to repeat itself

            energy_calc_simple(correct_output, a+1, falling_particles); //calculates "energy" for testing (falling atoms)

            printf("\n\tLast accepted energy = %lf\n\tCurrent energy = %lf", first_output->energy, correct_output->energy);
            num_for_energy_decision = num_gen(1, particle_list); /****/ //generates a random number between 0 and 1 for move acceptance/rejection
            ratio = exp((first_output->energy - correct_output->energy)/(boltzmann*temperature)); /****/ //calculates the value for comparison to num_for_energy_decision

            if(((first_output->energy) > (correct_output->energy)) || (ratio>*num_for_energy_decision))
            {
                copy_POSCAR(first_output, correct_output); //copying correct_output to first_output. correct_output is the most recent POSCAR that's had a random move successfully made. first_output is the POSCAR that gets written
                copy_POSCAR(second_output, correct_output); //copying correct_output to second_output. second_output is basically a holder POSCAR used for passing information
                write_POSCAR(first_output, POSCAR_file_path, replacement, 0);
                success = 1;
                accepted = 1;
                attempt++;
                printf("\n\n\t\t\tSUCCESS\n");
                write_XDATCAR(file_prefix, a_finished_on, particle_list, a+1); //appends the pre-existing XDATCAR file
            }
            else if((first_output->energy) < (correct_output->energy)) //if the new move is rejected, this SHOULD have the effect of undoing the move
            {
                a--;
                accepted = 0;
                attempt++;
                correct_output->rm_terminated = 1;
            }
            else if((first_output->energy)==(correct_output->energy)) //I separated this statement from the previous statement because of debugging paranoia, but it's the same
            {
                a--;
                attempt++;
                accepted = 0;
            }



            if(accepted==1)
            {
                if(correct_output->move_type==0)
                {
                    total_num_translate++;
                }
                else if(correct_output->move_type==1)
                {
                    total_num_rotate++;
                }

                if(a>=anneal_end)
                {
                    anneal = 0;
                }

                if((anneal==1) && (a>=anneal_begin)) //adjusts temperature linearly for simulated annealing option
                {
                    temperature-=((Tmax - Tmin)/(anneal_end-anneal_begin));
                }

                //printf("\nTemperature = %lf", temperature);
            }

            accepted_move_number+=accepted;
            VASP_energies(correct_output, attempt, correct_output->move_type, accepted, accepted_move_number, attempted_move_number, total_num_translate, total_num_rotate, 0, 0); //for use on Palmetto


            if(correct_output->rm_terminated==1) //This keeps track of how many times a random move has been terminated
            {
                printf("\n\tRandom_move() at POSCAR_%i has failed. \n\tIncrementing rm_termination_counter to %i and retrying the random move.", a+2, rm_termination_counter+1);
                rm_termination_counter++;
                //a--;
            }
        }

        if(rm_termination_counter>=rm_termination_limit) //if the current random move ATTEMPT has been terminated too much, this ends the loop in which random moves occur and the program wraps everything up
        {
            a_finished_on = a+1; //keeps track of the final posiiton of a (in counting numbers, not index numbers)
            a = max_iterations; //sets a to a value which breaks out of a loop
        }

        if(success==1) //If a move is successfully made, this resets the random move termination counter
        {
            rm_termination_counter = 0;
            success = 0;
        }

    } //this is the end of the loop that is broken out of if too many move attempts are made

    VASP_energies(correct_output, attempt, correct_output->move_type, accepted, accepted_move_number, attempted_move_number, total_num_translate, total_num_rotate, 1, trans_rotat_ratio);

    a_finished_on = a;

    if(rm_termination_counter>=rm_termination_limit)
    {
        printf("\n\n13) Random move procedure terminated due to excessive unsuccessful attempts\n\n");
    }


     //writes the final file, an XDATCAR, of the simulation
    printf("\n\n\nXDATCAR is written");
fclose(master_input);
    //write_POSCAR(first_output, POSCAR_file_path, particle_list, 1);
    printf("\nThe resubmittable POSCAR is written");
////////////////////////////////////////////////////////////////////

    //Frees memory and closes files

    free_placed_particles(particle_list);
    //free_placed_particles(replacement); //commented out because this fails if replacement is not properly filled due to an early termination

    printf("\n\n\n\n\t\t\tCREATE FREE_RESUBMIT_FILE()****************\n\n\n\n");
    free_POSCAR_struct(Initial);
    free_POSCAR_struct(first_output);
    free_POSCAR_struct(second_output);
    free_POSCAR_struct(correct_output);
    free(file_path);
    free(POSCAR_file_path);


    return 0;
}


  /////////////////////////////////
 ////    Read input files     ////
/////////////////////////////////


int read_poscar(POSCAR_struct* input, char* path, int resubmit) //must pre-allocate input
{
    char buffer[100];

    int a, b, c, d; //counter variables
    int num_types = input->num_types;
    int par_count = 0;

    char ident_holder[num_types][3]; //NumCols = 3 b/c two characters per symbol plus a terminating null character

    FILE* poscar = fopen(path,"r");

    if(poscar==NULL)
    {
        printf("\nCRASH: Failed to open POSCAR.\nMake sure all file names and file prefixes are correct.\n\n");
        exit(2);
    }


////////////////////////////////////////////////////////


    fgets(buffer, 100, poscar); //reads POSCAR header so fscanf picks up the lattice scale factor next

    fscanf(poscar, "%lf", &input->scale_factor);


    for(a=0; a<9; a++)
    {
        fscanf(poscar, "%lf", &input->latt_vec[a]);
    }



//loop the next line for num_types iterations
    for(d=0; d<num_types; d++)
    {
        fscanf(poscar, "%s", &input->par_ident[d]);
    }




    input->par_num = (int *) malloc(num_types*sizeof(int));

    for(b=0; b<num_types; b++)
    {
        fscanf(poscar, "%i", &input->par_num[b]);
        par_count+=input->par_num[b];
    }
    input->par_count = par_count;


////////////////////////////////////////////////////////////////

//Reads lines so fscanf can pick up the atom coordinates

    fgets(buffer, 100, poscar);

    fscanf(poscar,"\n%s", buffer); //Selective

    if((strcmp(buffer, "Selective")==0 || (strcmp(buffer, "selective"))==0))
    {
        fscanf(poscar, "%s", buffer); //dynamics
        //printf("\n3) %s\n", buffer);
    }

    if((strcmp(buffer, "dynamics")==0 || (strcmp(buffer, "Dynamics"))==0))
    {
        fscanf(poscar, "\n%s", buffer); //Direct
        //printf("\n4) %s\n", buffer);
    }

//////////////////////////////////////////////////////////////

    input->par_pos = (double*) malloc(3*par_count*sizeof(double)); //3 -> x, y, and z coordinates
    input->molecular_ID = (int*) malloc(par_count*sizeof(int));
    input->molecular_type_tag = (int*) malloc(par_count*sizeof(int));
    input->atom_ID = (int*) malloc(par_count*sizeof(int));
    input->bonded_to = (int*) malloc(2*par_count*sizeof(int));
    input->angle = (int*) malloc(3*par_count*sizeof(int));

    if(resubmit==0)
    {
        for(c=1; c<=par_count; c++)
        {
            fscanf(poscar, "%lf %lf %lf", &input->par_pos[3*c-3], &input->par_pos[3*c-2], &input->par_pos[3*c-1]);
            fgets(buffer, 100, poscar);
            input->molecular_ID[c-1] = -1;
        }
    }
    else if(resubmit==1)
    {
        for(c=1; c<=par_count; c++)
        {
            fscanf(poscar, "%lf %lf %lf #%i %i", &input->par_pos[3*c-3], &input->par_pos[3*c-2], &input->par_pos[3*c-1], &input->molecular_ID[c-1], &input->molecular_type_tag[c-1]);
            //printf("\nmolecular_ID = %i\tmolecular_type_tag = %i", input->molecular_ID[c-1], input->molecular_type_tag[c-1]);
            fgets(buffer, 100, poscar);
        }
    }


    fclose(poscar);

    return 0;
}


placed_particles* read_internal_coords(FILE* master_input, char* file_prefix)
{

    /*

        This function reads structure data for each kind of molecule into the molecule_coords and molecule_struct
        structures.
        There is one placed_particles structure (until replacement begins), as many instances of molecule_coords as there are kinds of molecules,
        and as many instances of molecule_struct as there are kinds of molecules.

    */

    placed_particles* molec_array = (placed_particles*) malloc(sizeof(placed_particles));

    FILE* fh = NULL;
    char* fh_name = (char*) malloc(50*sizeof(char));
    char* fh_filepath = NULL;
    int length = 0;
    int num_atoms = 0;
    int num_atom_types = 0;
    int* num_each_type = NULL;
    int* max_num_each_type = NULL;
    int max_number = 0;

    int a, b, num_types, d, fh_counter, e; //counter variables

    char* type_pointer = NULL;
    char* type_reset = NULL;

    int index_num_types = 0;

    char* ident_pointer = NULL;
    int* bonded_to = NULL;
    double* lengths_bond = NULL;
    int* third_atom = NULL;
    double* angles_bond = NULL;
    int* fourth_atom = NULL;
    double* dihedral = NULL;



///////////////////////////////////////////////////////////////////

    fscanf(master_input, "\nNumber of different types of fluid molecules: %i", &molec_array->num_types_molecules);
    //printf("\nmp_num_types = %i\n", molec_array->num_types_molecules);



    molec_array->molecule_info = (molecule_coords*) malloc(molec_array->num_types_molecules*sizeof(molecule_coords));

    //Begin reading information from pdb file (begin loop here)
    for(fh_counter = 0; fh_counter<molec_array->num_types_molecules; fh_counter++)
    {
        fscanf(master_input,"%i", &molec_array->molecule_info[fh_counter].max_number);

        molec_array->molecule_info[fh_counter].structure = (molecule_struct*) malloc(sizeof(molecule_struct));
        molec_array->molecule_info[fh_counter].replacing = 0;
        molec_array->molecule_info[fh_counter].molecular_type_tag = fh_counter;
        molec_array->terminate_replacing = 0;


        //read in pdb file name and open file

        fscanf(master_input,"%s", fh_name);

        length+=(strlen(fh_name) + strlen(file_prefix) + 1);

        fh_filepath = (char*) malloc(length*sizeof(char));

        sprintf(fh_filepath, "%s%s", file_prefix, fh_name);
printf("\nfh_filepath = %s", fh_filepath);
        fh = fopen(fh_filepath, "r");


        if(fh==NULL)
        {
            printf("\n\t\tFailed to open fh file!\n");
            exit(3);
        }

        //read in num_atoms

        fscanf(fh, "%i", &molec_array->molecule_info[fh_counter].structure->num_atoms);

        num_atoms = molec_array->molecule_info[fh_counter].structure->num_atoms;

//printf("\nnum_atoms = %i", molec_array->molecule_info[fh_counter].structure->num_atoms);


        //allocate memory for all pointers in molecule_struct and assign correct addresses to shortened pointers

        molec_array->molecule_info[fh_counter].structure->bonded_to = (int*) malloc(num_atoms*sizeof(int));
        molec_array->molecule_info[fh_counter].structure->lengths_bond = (double*) malloc(num_atoms*sizeof(double));
        molec_array->molecule_info[fh_counter].structure->third_atom = (int*) malloc(num_atoms*sizeof(int));
        molec_array->molecule_info[fh_counter].structure->angles_bond = (double*) malloc(num_atoms*sizeof(double));
        molec_array->molecule_info[fh_counter].structure->fourth_atom = (int*) malloc(num_atoms*sizeof(int));
        molec_array->molecule_info[fh_counter].structure->dihedral = (double*) malloc(num_atoms*sizeof(double));

        molec_array->molecule_info[fh_counter].structure->bond_angle_ranges = (double*) malloc(num_atoms*sizeof(double));

        for(a=0; a<num_atoms; a++)
        {
            molec_array->molecule_info[fh_counter].structure->bond_angle_ranges[a] = 0;
        }


        ident_pointer = molec_array->molecule_info[fh_counter].structure->atom_ident;
        bonded_to = molec_array->molecule_info[fh_counter].structure->bonded_to;
        lengths_bond = molec_array->molecule_info[fh_counter].structure->lengths_bond;
        third_atom = molec_array->molecule_info[fh_counter].structure->third_atom;
        angles_bond = molec_array->molecule_info[fh_counter].structure->angles_bond;
        fourth_atom = molec_array->molecule_info[fh_counter].structure->fourth_atom;
        dihedral = molec_array->molecule_info[fh_counter].structure->dihedral;

        //fill in members of molecule_struct
            //simultaneously read in atom_ident to molecule_struct, determine num_atom_types and num_each_type,
            //and fill in the rest of molecule_struct

        for(a=0; a<num_atoms; a++)
        {
            if(a==0)
            {
                fscanf(fh, "\n%s %i", molec_array->molecule_info[fh_counter].structure->atom_ident[a], &bonded_to[a]);

                //printf("\n\natom_ident[%i] = %s", a, molec_array->molecule_info[fh_counter].structure->atom_ident[a]);
                //printf("\n\tbonded_to[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->bonded_to[a]);
            }
            else if(a==1)
            {
                fscanf(fh, "\n%s %i %lf", molec_array->molecule_info[fh_counter].structure->atom_ident[a], &bonded_to[a], &lengths_bond[a]);

                //printf("\n\natom_ident[%i] = %s", a, molec_array->molecule_info[fh_counter].structure->atom_ident[a]);
                //printf("\n\tbonded_to[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->bonded_to[a]);
                //printf("\n\tlengths_bond[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->lengths_bond[a]);
            }
            else if(a==2)
            {
                fscanf(fh, "\n%s %i %lf %i %lf", molec_array->molecule_info[fh_counter].structure->atom_ident[a], &bonded_to[a], &lengths_bond[a], &third_atom[a], &angles_bond[a]);

                if(angles_bond[a]<0)
                {
                    fscanf(fh, "%lf", &molec_array->molecule_info[fh_counter].structure->bond_angle_ranges[a]);
                }


                //printf("\n\natom_ident[%i] = %s", a, molec_array->molecule_info[fh_counter].structure->atom_ident[a]);
                //printf("\n\tbonded_to[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->bonded_to[a]);
                //printf("\n\tlengths_bond[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->lengths_bond[a]);
                //printf("\n\tthird_atom[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->third_atom[a]);
                //printf("\n\tangles_bond[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->angles_bond[a]);
            }
            else if(a>=3)
            {
                fscanf(fh, "\n%s %i %lf %i %lf %i %lf", molec_array->molecule_info[fh_counter].structure->atom_ident[a], &bonded_to[a], &lengths_bond[a], &third_atom[a], &angles_bond[a], &fourth_atom[a], &dihedral[a]);

                //printf("\n\natom_ident[%i] = %s", a, molec_array->molecule_info[fh_counter].structure->atom_ident[a]);
                //printf("\n\tbonded_to[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->bonded_to[a]);
                //printf("\n\tlengths_bond[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->lengths_bond[a]);
                //printf("\n\tthird_atom[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->third_atom[a]);
                //printf("\n\tangles_bond[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->angles_bond[a]);
                //printf("\n\tfourth_atom[%i] = %i", a, molec_array->molecule_info[fh_counter].structure->fourth_atom[a]);
                //printf("\n\tdihedral[%i] = %lf", a, molec_array->molecule_info[fh_counter].structure->dihedral[a]);
            }
            else
            {
                printf("\nSomething went terribly, terribly wrong!\n");
                exit(4);
            }
        }
//printf("\n\n\n////////////////////////////////////////////////////\n\n\n");
        //Filling in type_ident[], num_atom_types, and num_each_type


        num_types = 1; //counter variable for type_ident array

        type_pointer = molec_array->molecule_info[fh_counter].type_ident;
        type_reset = molec_array->molecule_info[fh_counter].type_ident;

        ident_pointer = molec_array->molecule_info[fh_counter].structure->atom_ident;

        copyString(type_pointer, ident_pointer);
        //printf("\nMatch: type = %s\tident = %s",type_pointer, ident_pointer);
        ident_pointer+=3;

        molec_array->molecule_info[fh_counter].num_each_type = (int*) malloc(50*sizeof(int));
        molec_array->molecule_info[fh_counter].max_num_each_type = (int*) malloc(50*sizeof(int));

        max_number = molec_array->molecule_info[fh_counter].max_number;
        max_num_each_type = molec_array->molecule_info[fh_counter].max_num_each_type;
        num_each_type = molec_array->molecule_info[fh_counter].num_each_type;
        num_each_type[0] = 1;

        for(b=1; b<50; b++)
        {
            num_each_type[b] = 0;
        }



        for(b=0; b<(num_atoms-1); b++)
        {
            //printf("\n\n b = %i",b);


            for(d=0; d<num_types; d++)
            {
                if(strcmp(type_pointer, ident_pointer)!=0) //case if type and atom_ident aren't the same
                {
                    //printf("\nMismatch: type = %s\tident = %s", type_pointer, ident_pointer);
                    type_pointer+=3;
                    index_num_types++;

                    if(d==(num_types-1))
                    {
                        copyString(type_pointer, ident_pointer);
                        //printf("\tNew type = %s", type_pointer);
                        num_each_type[index_num_types]+=1;
                        num_types++;
                        d = num_types;
                    }
                }
                else if(strcmp(type_pointer, ident_pointer)==0)
                {
                    //printf("\nMatch: type = %s\tident = %s", type_pointer, ident_pointer);
                    d = num_types;
                    num_each_type[index_num_types]+=1;
                }
            }


            ident_pointer+=3;
            //printf("\n Advance ident_pointer\n");
            type_pointer = type_reset;
            index_num_types = 0;
        }

        molec_array->molecule_info[fh_counter].num_atom_types = num_types;


        //printf("\n\n\n");
        type_pointer = type_reset;

        for(b=0; b<num_types; b++)
        {
            max_num_each_type[b] = max_number*num_each_type[b];

            //printf("\ntype_pointer = %s\t max_num_each_type = %i", type_pointer, max_num_each_type[b]);
            //type_pointer+=3;
        }


        free(fh_filepath);
        fclose(fh);
    }

    free(fh_name);

    return molec_array;
}

  ////////////////////////////////////////////
 ////    Particle Placement Functions    ////
////////////////////////////////////////////


double* box_max(double sp_pos_carte[], int sp_count, POSCAR_struct* input)
{
    double* max = (double*) malloc(3*sizeof(double));

        if(max==NULL)
        {
            printf("\n\nFailed to allocate memory for max!\n\n");
            exit(17);
        }

    int a;
    int b;
    int n;

    double scale = 1.05; // keep this value greater than 1

    double l1 = (input->scale_factor)*(input->latt_vec[0]);
    double l2 = (input->scale_factor)*(input->latt_vec[1]);
    double l3 = (input->scale_factor)*(input->latt_vec[2]);
    double l4 = (input->scale_factor)*(input->latt_vec[3]);
    double l5 = (input->scale_factor)*(input->latt_vec[4]);
    double l6 = (input->scale_factor)*(input->latt_vec[5]);
    double l7 = (input->scale_factor)*(input->latt_vec[6]);
    double l8 = (input->scale_factor)*(input->latt_vec[7]);
    double l9 = (input->scale_factor)*(input->latt_vec[8]);

    max[0]=0;
    max[1]=0;
    max[2]=0;

    //defining x_max
    if(l1>0)  //if x lattice vector points in the positive x direction
    {
        if((l4<=0) && (l7<=0))
        {
            max[0] = scale*l1;
        }
        else if((l4>0) && (l7<=0))
        {
            max[0] = l4 + scale*l1;
        }
        else if((l7>0) && (l4<=0))
        {
            max[0] = l7 + scale*l1;
        }
        else if((l4>0) && (l7>0))
        {
            max[0] = l4 + l7 + scale*l1;
        }
    }
    else if(l1<0) //if x lattice vector points in the negative x direction
    {
        if((l4>=0) && (l7>=0))
        {
            max[0] = l4 + l7 - (scale-1)*l1;
        }
        else if((l4<0) && (l7>=0))
        {
            max[0] = l7 - (scale-1)*l1;
        }
        else if((l7<0) && (l4>=0))
        {
            max[0] = l4 - (scale-1)*l1;
        }
        else if((l4<0) && (l7<0))
        {
            max[0] = -(scale-1)*l1;
        }
    }

    //defining y_max
    if(l5>0) //if y lattice vector points in the positive y direction
    {
        if((l2<=0) && (l8<=0))
        {
            max[1] = scale*l5;
        }
        else if((l2>0) && (l8<=0))
        {
            max[1] = l2 + scale*l5;
        }
        else if((l8>0) && (l2<=0))
        {
            max[1] = l8 + scale*l5;
        }
        else if((l2>0) && (l8>0))
        {
            max[1] = l2 + l8 + scale*l5;
        }
    }
    else if(l5<0) //if y lattice vector points in the negative y direction
    {
        if((l2>=0) && (l8>=0))
        {
            max[1] = l2 + l8 - (scale-1)*l5;
        }
        else if((l2<0) && (l8>=0))
        {
            max[1] = l8 - (scale-1)*l5;
        }
        else if((l8<0) && (l2>=0))
        {
            max[1] = l2 - (scale-1)*l5;
        }
        else if((l2<0) && (l8<0))
        {
            max[1] = -(scale-1)*l5;
        }
    }

    //defining z_max
    if(l9>0) //if z lattice vector points in the positive z direction
    {
        if((l3<=0) && (l6<=0))
        {
            max[2] = scale*l9;
        }
        else if((l3>0) && (l6<=0))
        {
            max[2] = l3 + scale*l9;
        }
        else if((l6>0) && (l3<=0))
        {
            max[2] = l6 + scale*l9;
        }
        else if((l3>0) && (l6>0))
        {
            max[2] = l3 + l6 + scale*l9;
        }
    }
    else if(l9<0) //if z lattice vector points in the negative z direction
    {
        if((l3>=0) && (l6>=0))
        {
            max[2] = l3 + l6 - (scale-1)*l9;
        }
        else if((l3<0) && (l6>=0))
        {
            max[2] = l6 - (scale-1)*l9;
        }
        else if((l6<0) && (l3>=0))
        {
            max[2] = l3 - (scale-1)*l9;
        }
        else if((l3<0) && (l6<0))
        {
            max[2] = -(scale-1)*l9;
        }
    }

    //printf("\nbox_max in fxn: %f\t%f\t%f\n\n",max[0],max[1],max[2]);

    return max;
}


double* box_min(double sp_pos_carte[], int sp_count, POSCAR_struct* input)
{//For uneven surfaces, I will need to make the box_min values be 0 and rely on the placement checking routine for even coverage without gaps in the surface's valleys.
    double* min=(double*) malloc(3*sizeof(double));

    if(min==NULL)
    {
        printf("\n\nFailed to allocate memory for min!\n\n");
        exit(17);
    }

    //int a;
    int b;
    int n;
    //int c;

    double scale = 1.05;

    double l1 = input->scale_factor*input->latt_vec[0];
    double l2 = input->scale_factor*input->latt_vec[1];
    double l3 = input->scale_factor*input->latt_vec[2];
    double l4 = input->scale_factor*input->latt_vec[3];
    double l5 = input->scale_factor*input->latt_vec[4];
    double l6 = input->scale_factor*input->latt_vec[5];
    double l7 = input->scale_factor*input->latt_vec[6];
    double l8 = input->scale_factor*input->latt_vec[7];
    double l9 = input->scale_factor*input->latt_vec[8];

    //defining x_min
    if(l1>0) //if x lattice vector is pointing in the positive x direction
    {
        if((l4<=0) && (l7<=0))
        {
            min[0] = l4 + l7 - (scale-1)*l1;
        }
        else if((l4>0) && (l7<=0))
        {
            min[0] = l7 - (scale-1)*l1;
        }
        else if((l7>0) && (l4<=0))
        {
            min[0] = l4 - (scale-1)*l1;
        }
        else if((l7>0) && (l4>0))
        {
            min[0] = -(scale-1)*l1;
        }
    }
    else if(l1<0) //if x lattice vector is pointing in the negative x direction
    {
        if((l4>=0) && (l7>=0))
        {
            min[0] = scale*l1;
        }
        else if((l4<0) && (l7>=0))
        {
            min[0] = l4 + scale*l1;
        }
        else if((l7>0) && (l4<=0))
        {
            min[0] = l7 + scale*l1;
        }
        else if((l7>0) && (l4>0))
        {
            min[0] = l7 + l4 + scale*l1;
        }
    }


    //defining y_min
    if(l5>0) //if y lattice vector is pointing in the positive y direction
    {
        if((l2<=0) && (l8<=0))
        {
            min[1] = l2 + l8 - (scale-1)*l5;
        }
        else if((l2>0) && (l8<=0))
        {
            min[1] = l8 - (scale-1)*l5;
        }
        else if((l8>0) && (l2<=0))
        {
            min[1] = l2 - (scale-1)*l5;
        }
        else if((l8>0) && (l2>0))
        {
            min[1] = -(scale-1)*l5;
        }
    }
    else if(l5<0) //if y lattice vector is pointing in the negative y direction
    {
        if((l2>=0) && (l8>=0))
        {
            min[1] = scale*l5;
        }
        else if((l2<0) && (l8>=0))
        {
            min[1] = l2 + scale*l5;
        }
        else if((l8>0) && (l2<=0))
        {
            min[1] = l8 + scale*l5;
        }
        else if((l8>0) && (l2>0))
        {
            min[1] = l8 + l2 + scale*l5;
        }
    }

    //defining z_min
    if(l9>0) //if z lattice vector is pointing in the positive z direction
    {
        if((l3<=0) && (l6<=0))
        {
            min[2] = l3 + l6 - (scale-1)*l9;
        }
        else if((l3>0) && (l6<=0))
        {
            min[2] = l6 - (scale-1)*l9;
        }
        else if((l6>0) && (l3<=0))
        {
            min[2] = l3 - (scale-1)*l9;
        }
        else if((l6>0) && (l3>0))
        {
            min[2] = -(scale-1)*l9;
        }
    }
    else if(l9<0) //if z lattice vector is pointing in the negative z direction
    {
        if((l3>=0) && (l6>=0))
        {
            min[2] = scale*l9;
        }
        else if((l3<0) && (l6>=0))
        {
            min[2] = l3 + scale*l9;
        }
        else if((l6>0) && (l3<=0))
        {
            min[2] = l6 + scale*l9;
        }
        else if((l6>0) && (l3>0))
        {
            min[2] = l6 + l3 + scale*l9;
        }
    }

    //printf("\nbox_min in fxn: %lf\t%lf\t%lf",min[0],min[1],min[2]);

    return min;
}


int particle_placement_shell(placed_particles* particle_list, POSCAR_struct* input, int zone, double min_scale, double max_scale, double z_min, double z_max, double maximum_translation)
{
    /**

    -Call box-sizing function(s)
    -Place particle atom by atom
        -Check placement of each atom
            -If placement fails, fall back to first atom with a degree of freedom and replace it
            -If this fails a specified number of times, move back to a previously placed atom with a degree of freedom
            -If the program falls all the way back to the first atom and still fails, begin placing the next type of particle

    **/

    int rm_terminated = 0;

    double* sp_pos_carte = frac_to_cartesian(input->latt_vec, input->par_pos, input->par_count, input->scale_factor);
    double* max_box_dim = box_max(sp_pos_carte, input->par_count, input);
    double* min_box_dim = box_min(sp_pos_carte, input->par_count, input);


    double* coords_holder = NULL;
    double* coords_pointer = NULL;

    int zz;
    double max_dimension = 0;

    for(zz=0; zz<3; zz++)
    {
        if(max_box_dim[zz]>max_dimension)
        {
            max_dimension = max_box_dim[zz];
        }
    }


    /**Information, counters, and pointers needed to place each kind of molecule**/
    int terminate = 0; //if set to 1, terminates placement of a type of molecule if certain criteria are met
                            //see next three lines of grey text for criteria

    /**IMPORTANT NOTE--> in placement_check function, an atom's degrees of freedom must be taken into account
                        -If an atom's information contains a FIXED dihedral angle and its placement has failed,
                            the placement routine must retreat to the 3rd atom (first atom w/o a dihedral angle)**/



    int num_types_molecules = particle_list->num_types_molecules; //upper limit for current_molec counter variable
    int current_molec; //counter variable for outer loop

    int max_number_molecules = 0;

    int max_number_atoms = 0; //maximum allowed number of the current molecule type
                    //This will be the upper limit to the middle loop
    int actual_number_placed = 0; //Records the actual number of the current molecule type placed
                                  //counter variable for the middle loop

    int num_atoms = 0; //The number of atoms in every molecule of the current molecule type
                       //This will be the upper limit to a counter variable in the innermost loop
    int current_atom = 1; //counter variable

    int attempt = 0; //number of attempted placements for the current atom.
                    /**The upper limit is currently    **/

    int nth_atom_reset_counter = 0;
    int third_atom_reset_counter = 0;
    int second_atom_reset_counter = 0;
    int first_atom_reset_counter = 0;

    double success = 0;

    int one = 1; //don't delete these! Maybe? I'm using these because I'm not sure if integers that aren't a defined variable are declared as double/float or int by the compiler
    int zero = 0;


    molecule_struct* structure_info = NULL;

    int* bonded_to = NULL;
    double* lengths_bond = NULL;
    int* third_atom = NULL;
    double* angles_bond = NULL;
    int* fourth_atom = NULL;
    double* dihedral = NULL;

    int molecule_success = 0;

    particle_list->molecule_count = 0;

    int molec_ID_index;

    if(zone==1)
    {
        min_box_dim[2] = z_min;
    }
    else if(zone==2)
    {
        max_box_dim[2] = z_max;
    }
    else if(zone==3)
    {
        min_box_dim[2] = z_min;
        max_box_dim[2] = z_max;
    }

    for(current_molec=0; current_molec<num_types_molecules; current_molec++)
    {
        particle_list->molecule_info[current_molec].actual_number_placed = 0;
    }

    int a, aa, b, c, d, e;
    int rotation_failure = 0;
    int translation_failure = 0;

    /// Translation coordinate arrays
    double* t_original_pos = (double*) malloc(3*num_atoms*sizeof(double));
    double* t_original_pos_cartesian = (double*) malloc(3*num_atoms*sizeof(double));
    double* t_translated_pos = (double*) malloc(3*num_atoms*sizeof(double));


    int current_atom_ID = 0;
    int first_atom_ID = 0;


    for(current_molec=0; current_molec<num_types_molecules; current_molec++) //loops through all molecule types
    {
        particle_list->molecule_info[current_molec].actual_number_placed = 0;

        /**Go through all types of molecules here**/
        max_number_molecules = particle_list->molecule_info[current_molec].max_number;
        max_number_atoms = (max_number_molecules)*(particle_list->molecule_info[current_molec].structure->num_atoms);
        structure_info = particle_list->molecule_info[current_molec].structure;
        particle_list->molecule_info[current_molec].coords = (double*) malloc(3*max_number_atoms*sizeof(double));
        particle_list->molecule_info[current_molec].molecular_ID = (int*) malloc(max_number_atoms*sizeof(double));
        coords_pointer = particle_list->molecule_info[current_molec].coords;

        particle_list->molecule_info[current_molec].atom_ID = (int*) malloc(max_number_atoms*sizeof(int));
        particle_list->molecule_info[current_molec].bonded_to = (int*) malloc(2*max_number_atoms*sizeof(int));
        particle_list->molecule_info[current_molec].angle = (int*) malloc(3*max_number_atoms*sizeof(int));

        particle_list->molecule_info[current_molec].min_distance_array = (double*) malloc(particle_list->molecule_info[current_molec].structure->num_atoms*sizeof(double));



        while((actual_number_placed<max_number_molecules) && (terminate==0) ) //loops through all instances of a molecule type
        { /**Going through all molecules of a particular type**/



        bonded_to = structure_info->bonded_to;
        lengths_bond = structure_info->lengths_bond;
        third_atom = structure_info->third_atom;
        angles_bond = structure_info->angles_bond;
        fourth_atom = structure_info->fourth_atom;
        dihedral = structure_info->dihedral;

        num_atoms = structure_info->num_atoms;

        double position_and_info[4] = {0};

    /// Rotation coordinate arrays
        double* r_original_pos = (double*) malloc(3*num_atoms*sizeof(double)); //original atom positions (fractional coordinates)
        double* r_original_pos_cartesian; //original atom positions (cartesian coordinates)
        double* r_new_pos = (double*) malloc(3*num_atoms*sizeof(double)); //newly generated atom positions
        double* r_rotated_pos = (double*) malloc(3*num_atoms*sizeof(double)); //correctly rotated atom positions (correctly rotated doesn't mean it's been checked yet. It's just a rotation)
        double* r_translation_vector = (double*) malloc(3*sizeof(double)); //a vector pointing from an old atom position to a new atom position in the rotated molecule
        int rotate_about;


        while((current_atom<=num_atoms) && (terminate==0) && (molecule_success==0)) //loops through all atoms in a molecule type
        { /**Placing all atoms in a particular molecule**/
            translation_failure = 0;
            rotation_failure = 0;

            while((terminate==0) && (isEqual(zero, success))) //loop focusing on placement of current atom in current molecule type
            {
                coords_holder = atom_placer(input,
                                            coords_pointer,
                                            particle_list,
                                            &particle_list->molecule_info[current_molec],
                                            current_atom,
                                            bonded_to[current_atom-1],
                                            lengths_bond[current_atom-1],
                                            third_atom[current_atom-1],
                                            angles_bond[current_atom-1],
                                            fourth_atom[current_atom-1],
                                            dihedral[current_atom-1],
                                            max_box_dim,
                                            min_box_dim,
                                            min_scale,
                                            max_scale,
                                            z_min,
                                            z_max);
                success = coords_holder[0];


                if(((particle_list->molecule_info->rotate==1) || (particle_list->molecule_info->translate==1)) && (success==1))
                { //printf("\n\nreplaced_par_pos_indices[%i] = %i", current_atom-1, particle_list->molecule_info->replaced_par_pos_indices[current_atom-1]);
                    r_original_pos[3*current_atom-3] = input->par_pos[3*particle_list->molecule_info->replaced_par_pos_indices[current_atom-1]];
                    r_original_pos[3*current_atom-2] = input->par_pos[3*particle_list->molecule_info->replaced_par_pos_indices[current_atom-1]+1];
                    r_original_pos[3*current_atom-1] = input->par_pos[3*particle_list->molecule_info->replaced_par_pos_indices[current_atom-1]+2];

                    //printf("\nr_original_pos[%i] = %lf", 3*current_atom-3, r_original_pos[3*current_atom-3]);
                    //printf("\nr_original_pos[%i] = %lf", 3*current_atom-2, r_original_pos[3*current_atom-2]);
                    //printf("\nr_original_pos[%i] = %lf", 3*current_atom-1, r_original_pos[3*current_atom-1]);

                    r_new_pos[3*current_atom-3] = coords_holder[1];
                    r_new_pos[3*current_atom-2] = coords_holder[2];
                    r_new_pos[3*current_atom-1] = coords_holder[3];

                    //printf("\ncoords_holder[0] = %lf", coords_holder[1]);
                    //printf("\ncoords_holder[1] = %lf", coords_holder[2]);
                    //printf("\ncoords_holder[2] = %lf", coords_holder[3]);

                    //printf("\n\nr_new_pos[%i] = %lf", 3*current_atom-3, r_new_pos[3*current_atom-3]);
                    //printf("\nr_new_pos[%i] = %lf", 3*current_atom-2, r_new_pos[3*current_atom-2]);
                    //printf("\nr_new_pos[%i] = %lf", 3*current_atom-1, r_new_pos[3*current_atom-1]);
                }
                if(isEqual(one, success))
                {
                    coords_pointer[(3*((actual_number_placed*num_atoms)+current_atom))-3] = coords_holder[1]; //coords_pointer has to be filled so that the placement of subsequent new atoms in the placement routine can happen. otherwise, the new positions are all someting like -1.#IND0 (?)
                    coords_pointer[(3*((actual_number_placed*num_atoms)+current_atom))-2] = coords_holder[2];
                    coords_pointer[(3*((actual_number_placed*num_atoms)+current_atom))-1] = coords_holder[3];

                    //if(rotation_failure==0) // potential problem section, but it's also essential to keeping the molecular ID tags correct. Maybe it's not a potential problem section?
                    //{
                        molec_ID_index = ((particle_list->molecule_info[current_molec].actual_number_placed)*(particle_list->molecule_info[current_molec].structure->num_atoms))+current_atom-1;

                        particle_list->molecule_info[current_molec].molecular_ID[molec_ID_index] = particle_list->molecule_count + actual_number_placed;

                        if(particle_list->molecule_info[current_molec].replacing==1)
                        {
                            particle_list->molecule_info[current_molec].molecular_ID[molec_ID_index] = particle_list->molecule_info[current_molec].replaced_molecular_ID;
                        }
                    //}
                }

                if((isEqual(zero, success) && (current_atom>=4)))
                {
                    nth_atom_reset_counter++;
                    if(nth_atom_reset_counter>=nth_atom_limit)
                    {
                        current_atom = 3;
                        nth_atom_reset_counter = 0;
                        success = 1; //not actually succesful. Using this to bypass the following 'if' statement(s)
                    }
                }

                if(isEqual(zero, success) && (current_atom==3))
                {
                    third_atom_reset_counter++;
                    if(third_atom_reset_counter>=third_atom_limit)
                    {
                        current_atom = 2;
                        third_atom_reset_counter = 0;
                        success = 1; //not actually succesful. Using this to bypass the following 'if' statement(s)
                    }
                }

                if((isEqual(zero, success)) && (current_atom==2))
                {
                    second_atom_reset_counter++;
                    if(second_atom_reset_counter>=second_atom_limit)
                    {
                        first_atom_reset_counter++;
                        current_atom = 1;
                        second_atom_reset_counter = 0;
                        success = 1; //not actually succesful. Using this to bypass the following 'if' statement(s)
                    }
                }

                if(isEqual(zero, success) && (current_atom==1))
                {
                    first_atom_reset_counter++;
                    second_atom_reset_counter = 0;
                    if(first_atom_reset_counter>=first_atom_limit)
                    {
                        terminate = 1;
                        if(particle_list->molecule_info[0].replacing==1)
                        {
                            particle_list->terminate_replacing = 1;
                        }
                        printf("\nTerminating placement routine (replacing = %i)\n", particle_list->molecule_info[current_molec].replacing);
                    }
                }
            } ///end of loop

                if((current_atom==(num_atoms)) && (isEqual(one, success)))
                {
                    first_atom_reset_counter = 0;
                    molecule_success = 1;
                    //printf("\nmolecule_success = %i", molecule_success);
                }

                success = 0;
                if(isEqual(one, coords_holder[0]))
                {
                    if(current_atom==(num_atoms))
                    {
                        success = 1;
                    }
                    else if(current_atom<(num_atoms))
                    {
                        current_atom++;
                    }
                }
            free(coords_holder);
        }///end of loop

            if((particle_list->molecule_info->rotate==1) || (particle_list->molecule_info->translate==1))
            {
                r_original_pos_cartesian = frac_to_cartesian(input->latt_vec, r_original_pos, num_atoms, input->scale_factor);

                rotate_about = particle_list->molecule_info->rotate_about_this_atom;

                if(particle_list->molecule_info->rotate==1)
                {
                    r_translation_vector[0] = r_original_pos_cartesian[3*rotate_about] - r_new_pos[3*rotate_about];
                    r_translation_vector[1] = r_original_pos_cartesian[3*rotate_about+1] - r_new_pos[3*rotate_about+1];
                    r_translation_vector[2] = r_original_pos_cartesian[3*rotate_about+2] - r_new_pos[3*rotate_about+2];
                }
                else if(particle_list->molecule_info->translate==1)
                {
                    double* source = num_gen(4, particle_list);
                    double scale = source[0]*maximum_translation;

                    r_translation_vector[0] = (source[1]*2)-1;
                    r_translation_vector[1] = (source[2]*2)-1;
                    r_translation_vector[2] = (source[3]*2)-1;

                    unitize(r_translation_vector);

                    r_translation_vector[0] = r_translation_vector[0]*scale;
                    r_translation_vector[1] = r_translation_vector[1]*scale;
                    r_translation_vector[2] = r_translation_vector[2]*scale;
                }

                if(particle_list->molecule_info->rotate==1)
                {
                    for(b=0; b<num_atoms; b++)
                    {
                        r_rotated_pos[3*b] = r_new_pos[3*b] + r_translation_vector[0];
                        r_rotated_pos[3*b+1] = r_new_pos[3*b+1] + r_translation_vector[1];
                        r_rotated_pos[3*b+2] = r_new_pos[3*b+2] + r_translation_vector[2];
                    }
                }
                else if(particle_list->molecule_info->translate==1)
                {
                    for(b=0; b<num_atoms; b++)
                    {
                        r_rotated_pos[3*b] = r_original_pos_cartesian[3*b] + r_translation_vector[0];
                        r_rotated_pos[3*b+1] = r_original_pos_cartesian[3*b+1] + r_translation_vector[1];
                        r_rotated_pos[3*b+2] = r_original_pos_cartesian[3*b+2] + r_translation_vector[2];
                    }
                }


                b=0;
                molecule_success=1;
                while((b<num_atoms) && (molecule_success==1))
                {
                    current_atom = b+1;

                    position_and_info[1] = r_rotated_pos[3*b]; //This has to be filled so that the placement checking functions can work
                    position_and_info[2] = r_rotated_pos[3*b+1];
                    position_and_info[3] = r_rotated_pos[3*b+2];

                    molecule_success = placement_checker_POSCAR(input, particle_list->molecule_info, current_atom, min_scale, max_scale, position_and_info, z_min, z_max);

                    b++;
                }

                if(molecule_success==0)
                {
                    rm_terminated = 1;
                }
                else if(molecule_success==1)
                {
                    rm_terminated = 0;
                    for(c=0; c<num_atoms; c++)
                    {
                        coords_pointer[(3*((actual_number_placed*num_atoms)+c+1))-3] = r_rotated_pos[3*c];
                        coords_pointer[(3*((actual_number_placed*num_atoms)+c+1))-2] = r_rotated_pos[3*c+1];
                        coords_pointer[(3*((actual_number_placed*num_atoms)+c+1))-1] = r_rotated_pos[3*c+2];
                    }
                }
            }

            success = 0;
            actual_number_placed+=molecule_success;
            //printf("\nactual_number_placed = %i", actual_number_placed);
            particle_list->molecule_info[current_molec].actual_number_placed = actual_number_placed;
            molecule_success = 0;
            current_atom = 1;
        }///end of loop

        particle_list->molecule_count+=actual_number_placed;
        particle_list->molecule_info[current_molec].actual_number_placed = actual_number_placed;

        first_atom_ID = current_atom_ID;
        e = 0;

        for(c=0; c<actual_number_placed; c++)
        {
            for(d=0; d<particle_list->molecule_info[current_molec].structure->num_atoms; d++)
            {
                particle_list->molecule_info[current_molec].bonded_to[2*e] = first_atom_ID + d;
                particle_list->molecule_info[current_molec].bonded_to[2*e+1] = first_atom_ID + bonded_to[d]-1; //the -1 accounts for the bonded_to/third_atom arrays starting to count at 1 instead of 0

                if((d>1) && (particle_list->molecule_info[current_molec].structure->num_atoms>2))
                {
                    particle_list->molecule_info[current_molec].angle[3*e] = first_atom_ID + d;
                    particle_list->molecule_info[current_molec].angle[3*e+1] = first_atom_ID + bonded_to[d]-1; //the -1 accounts for the bonded_to/third_atom arrays starting to count at 1 instead of 0
                    particle_list->molecule_info[current_molec].angle[3*e+2] = first_atom_ID + third_atom[d]-1; //the -1 accounts for the bonded_to/third_atom arrays starting to count at 1 instead of 0
                }
                else
                {
                    particle_list->molecule_info[current_molec].angle[3*e] = -1;
                    particle_list->molecule_info[current_molec].angle[3*e+1] = -1;
                    particle_list->molecule_info[current_molec].angle[3*e+2] = -1;
                }

                e++;

            }

            first_atom_ID += particle_list->molecule_info[current_molec].structure->num_atoms;
        }

        for(c=0; c<(actual_number_placed*particle_list->molecule_info[current_molec].structure->num_atoms); c++)
        {
            particle_list->molecule_info[current_molec].atom_ID[c] = current_atom_ID;
            current_atom_ID++;
        }

        actual_number_placed = 0;
        terminate = 0;
    } ///end of loop

    free(sp_pos_carte);
    free(max_box_dim);
    free(min_box_dim);

    return rm_terminated;
}


double* atom_placer(POSCAR_struct* input, double* coords_pointer, placed_particles* par_list, molecule_coords* molec_info, int atom_number, int bonded_to, double lengths_bond, int third_atom, double angles_bond, int fourth_atom, double dihedral, double max_dim[], double min_dim[], double min_scale, double max_scale, int z_min, int z_max)
{// note: doesn't actually place atoms
    double* position_and_info = (double*) malloc(4*sizeof(double)); //{success?, x, y, y} , where success == (0 or 1), 0 = unsuccessful, 1= successful
    double* pos = NULL;

    int a,b;

    int num_atoms_placed = 0;

    for(b=0; b<molec_info->num_atom_types; b++)
    {
        num_atoms_placed += molec_info->num_each_type[b];
    }

    num_atoms_placed = num_atoms_placed*molec_info->actual_number_placed;

    if(atom_number==1)
    {
        //printf("\n--------------------Placing %i--------------------", atom_number);

        pos = place_atom(par_list, max_dim, min_dim);

        for(a=0; a<3; a++)
        {
            position_and_info[a+1] = pos[a];
            //printf("\npos1[%i] = %lf", a, pos[a]);
        }
    }
    else if(atom_number==2)
    {
        //generate a random vector, centered on the first atom, with length = lengths_bond[i], pointing to the 2nd atom's position

        //printf("\n--------------------Placing %i--------------------", atom_number);

        /** Call place_with_bond_length(). Needs: coords[], bonded_to, bond length**/

        pos = place_with_bond_length(par_list, coords_pointer, bonded_to+num_atoms_placed, lengths_bond);

        for(a=0; a<3; a++)
        {
            position_and_info[a+1] = pos[a];
            //printf("\npos2[%i] = %lf", a, pos[a]);
        }

    }
    else if(atom_number==3)
    {
        //places an atom randomly in a circle in 3D space, the equation of which is determined by a bond length, bond angle, and two atom positions

        //printf("\n--------------------Placing %i--------------------", atom_number);



        /** Call place_with_bond_angle(). Needs: coords[], bonded_to, bond length, third atom, bond angle**/

        pos = place_with_bond_angle(par_list, coords_pointer, bonded_to+num_atoms_placed, lengths_bond, third_atom+num_atoms_placed, angles_bond);

        for(a=0; a<3; a++)
        {
            position_and_info[a+1] = pos[a];
            //printf("\npos3[%i] = %lf", a, pos[a]);
        }
    }
    else if(atom_number>=4)
    {


        //places an atom using the previous three atom positions, a bond angle, a bond length, and a dihedral angle

        //printf("\n\n\n--------------------Placing %i--------------------", atom_number);
        //printf("\n\n\nWARNING: DIHEDRAL PLACEMENT ROUTINE NOT IMPLEMENTED\n\n\n");
        /** Call place_with_dihedral(). Needs: coords[], bonded_to, bond length, third atom, bond angle, fourth atom, dihedral**/

        pos = place_with_dihedral(par_list, coords_pointer, bonded_to+num_atoms_placed, lengths_bond, third_atom+num_atoms_placed, angles_bond, fourth_atom+num_atoms_placed, dihedral);

        for(a=0; a<3; a++)
        {
            position_and_info[a+1] = pos[a];
            //printf("\npos[%i] = %lf", a, pos[a]);
        }

    }
    else
    {
        printf("\nYa sent a screwy atom_number value ( %i ), ya moran. (The accepted range is [1, pos. inf.)\n", atom_number);
        exit(10);
    }

    /**Call placement checkers here to determine success**/

    if((molec_info->rotate==0) || (molec_info->replacing!=1))
    {
        position_and_info[0] = placement_checker_POSCAR(input, molec_info, atom_number, min_scale, max_scale, position_and_info, z_min, z_max);

        if(position_and_info[0]==1)
        {
            position_and_info[0] = placement_checker_placedparticles(input, par_list, molec_info, atom_number, min_scale, max_scale, position_and_info);
        }
    }
    else if(molec_info->rotate==1)
    {
        position_and_info[0] = 1;
    }

//printf("\natom_placer-- min_scale = %lf\tmax_scale = %lf", min_scale, max_scale);
    //printf("\n\n\tpos[x] = %lf\n\tpos[y] = %lf\n\tpos[z] = %lf", position_and_info[1], position_and_info[2], position_and_info[3]);
//printf("\nposition_and_info[0] = %lf", position_and_info[0]);

    free(pos);
    return position_and_info;
}


double* place_atom(placed_particles* par_list, double max_dim[], double min_dim[])
{
    double* pos = num_gen(3, par_list);
    int a;

    for(a=0; a<3; a++)
    {
        pos[a] =(pos[a]*(max_dim[a]-min_dim[a]))+min_dim[a];
        //printf("\nmax_dim[%i] = %lf\tmin_dim[%i] = %lf", a, max_dim[a], a, min_dim[a]);
    }

    //printf("\n\n\tx = %lf\n\ty = %lf\n\tz = %lf", pos[0], pos[1], pos[2]);


    //printf("\n--------------------PLACE ATOM--------------------");
    //printf("\nbonded_to = %i", (a+1)/3);


    return pos;
}


double* place_with_bond_length(placed_particles* par_list, double coords[], int bonded_to, double bond_length)
{
    double* pos = (double*) malloc(3*sizeof(double));
    double theta, phi; //must be in radians for trig functions to work properly
    double w, x, y, z;
    double* angles = num_gen(2, par_list);

    int a;

    theta = angles[0]*2*pi;
    phi = angles[1]*2*pi;

    z = bond_length*sin(phi);
    w = pow(((bond_length*bond_length) - (z*z)), 0.5);

    x = w*cos(theta);
    y = w*sin(theta);

    pos[0] = coords[(3*bonded_to)-3] + x;
    pos[1] = coords[(3*bonded_to)-2] + y;
    pos[2] = coords[(3*bonded_to)-1] + z;

    //printf("\n--------------------BOND LENGTH--------------------");
    //printf("\nbonded_to = %i", bonded_to);

    free(angles);

    return pos;
}


double* place_with_bond_angle(placed_particles* par_list, double coords[], int bonded_to, double bond_length, int third_atom, double bond_angle)
{
    double atom_distance, atom_angle;

    double* c_pos = (double*) malloc(3*sizeof(double));

    double* xy_holder = num_gen(2, par_list);

        xy_holder[0] = (xy_holder[0]*30)-15;
        xy_holder[1] = (xy_holder[1]*30)-15;

    double x,y,z;



//printf("\n--------------------BOND ANGLE--------------------");
//printf("\nthird_atom = %i\n\n", third_atom);


//This section calculates the vertex's coordinates
    double tail_pos[3] = {coords[(3*third_atom)-3], coords[(3*third_atom)-2], coords[(3*third_atom)-1]}; //b_pos
    double head_pos[3] = {coords[(3*bonded_to)-3], coords[(3*bonded_to)-2], coords[(3*bonded_to)-1]};  //a_pos

    double vertex[3] = {0,0,0};

    double* normal_vector = create_unit_vector(head_pos, tail_pos); //creates a unit vector parallel to the bond between atoms 1 and 2.

    double vertex_distance = bond_length*sin((bond_angle-90)*(pi/180)); //distance from second atom to vertex
    double vertex_to_third_atom = bond_length*cos((bond_angle-90)*(pi/180)); //distance from vertex to third atom

    vertex[0] = (normal_vector[0]*vertex_distance) + head_pos[0];
    vertex[1] = (normal_vector[1]*vertex_distance) + head_pos[1];
    vertex[2] = (normal_vector[2]*vertex_distance) + head_pos[2];


    double z_new = ((normal_vector[0]/normal_vector[2])*(vertex[0]-xy_holder[0])) + ((normal_vector[1]/normal_vector[2])*(vertex[1]-xy_holder[1])) + vertex[2];

    double new_point[3] = {xy_holder[0], xy_holder[1], z_new};

    double* vertex_unit = create_unit_vector(new_point, vertex);

    x = vertex[0] + (vertex_to_third_atom*vertex_unit[0]);
    y = vertex[1] + (vertex_to_third_atom*vertex_unit[1]);
    z = vertex[2] + (vertex_to_third_atom*vertex_unit[2]);


    c_pos[0] = x;
    c_pos[1] = y;
    c_pos[2] = z;


//Results testing

    double* bc_unit = create_unit_vector(c_pos, head_pos);
    double* ba_unit = create_unit_vector(tail_pos, head_pos);

    atom_angle = (180/(pi))*acos(dot(bc_unit, ba_unit));

    atom_distance = dist(head_pos[0], head_pos[1], head_pos[2], c_pos[0], c_pos[1], c_pos[2]);

    //printf("\natom_angle = %lf", atom_angle);
    //printf("\natom_distance = %lf", atom_distance);


    //variables for checking ab_angle
/*
    double* ab_unit_vector = create_unit_vector(tail_pos, head_pos);

    double ab_angle = atan(ab_unit_vector[1]/ab_unit_vector[0]) * (180/pi);

    if ((ab_unit_vector[0]<0 && ab_unit_vector[1]>=0) || (ab_unit_vector[0]<0 && ab_unit_vector[1]<0))  //adjusts for limited range of atan() function
    {
        //printf("\nbc_angle = %lf", bc_angle);
        ab_angle+=180;
        printf("\nab_angle = %lf", ab_angle);
    }

    if(ab_angle<0)
    {
        ab_angle+=360;
    }

    printf("\nab_angle = %lf", ab_angle);

    printf("\n\na_pos[0] = %lf", head_pos[0]);
    printf("\na_pos[1] = %lf", head_pos[1]);
    printf("\na_pos[2] = %lf", head_pos[2]);

    printf("\n\nb_pos[0] = %lf", tail_pos[0]);
    printf("\nb_pos[1] = %lf", tail_pos[1]);
    printf("\nb_pos[2] = %lf", tail_pos[2]);
*/

    free(bc_unit);
    free(ba_unit);
    free(vertex_unit);
    free(normal_vector);
    free(xy_holder);

    return c_pos;
}


double* place_with_dihedral(placed_particles* par_list, double coords[], int bonded_to, double bond_length, int third_atom, double bond_angle, int fourth_atom, double dihedral)
{

    ///rules: n1.n2 = cos(dihedral), n2.b = 0, ||n2|| = 1

    double point_three[3] = {coords[3*fourth_atom-3], coords[3*fourth_atom-2], coords[3*fourth_atom-1]};
    double point_two[3] = {coords[3*third_atom-3], coords[3*third_atom-2], coords[3*third_atom-1]};
    double point_one[3] = {coords[3*bonded_to-3], coords[3*bonded_to-2], coords[3*bonded_to-1]};
    double* point_four = (double*) malloc(3*sizeof(double)); //this is the unknown position

    double vector_a[3] = {point_three[0]-point_two[0], point_three[1]-point_two[1], point_three[2]-point_two[2]};//{point_two[0]-point_three[0], point_two[1]-point_three[1], point_two[2]-point_three[2]};
    double vector_b[3] = {point_two[0]-point_one[0], point_two[1]-point_one[1], point_two[2]-point_one[2]};

    double length_b = dist(point_one[0], point_one[1], point_one[2], point_two[0], point_two[1], point_two[2]);

    double zerotol = 0.00001; //Is this a good value? Maybe it should be smaller. Maybe it doesn't matter.

    double normal_vector_two[3] = {0}; //this is the unknown normal vector
    double* normal_vector_one = cross(vector_b, vector_a);
    unitize(normal_vector_one);

    double degtorad = pi/180;

    double thetaabc = dihedral*degtorad;
    double costhetaabc = cos(thetaabc);
    double neg = 1;

    double n22, n21, n20;

    double sqrtterm1, sqrtterm2;

    if(thetaabc <= pi)
    {
        neg = 1;
    }
    else
    {
        neg = -1;
    }

    sqrtterm1 = -(vector_b[0]*vector_b[0])*(costhetaabc*costhetaabc) + (vector_b[0]*vector_b[0])*(normal_vector_one[1]*normal_vector_one[1]) + (vector_b[0]*vector_b[0])*(normal_vector_one[2]*normal_vector_one[2]) - 2*vector_b[0]*vector_b[1]*normal_vector_one[0]*normal_vector_one[1] - 2*vector_b[0]*vector_b[2]*normal_vector_one[0]*normal_vector_one[2] - vector_b[1]*vector_b[1]*costhetaabc*costhetaabc + vector_b[1]*vector_b[1]*normal_vector_one[0]*normal_vector_one[0] + vector_b[1]*vector_b[1]*normal_vector_one[2]*normal_vector_one[2] - 2*vector_b[1]*vector_b[2]*normal_vector_one[1]*normal_vector_one[2] - vector_b[2]*vector_b[2]*costhetaabc*costhetaabc + vector_b[2]*vector_b[2]*normal_vector_one[0]*normal_vector_one[0] + vector_b[2]*vector_b[2]*normal_vector_one[1]*normal_vector_one[1];
    sqrtterm2 = -(vector_b[0]*vector_b[0])*(costhetaabc*costhetaabc) + (vector_b[0]*vector_b[0])*(normal_vector_one[1]*normal_vector_one[1]) + (vector_b[0]*vector_b[0])*(normal_vector_one[2]*normal_vector_one[2]) - 2*vector_b[0]*vector_b[1]*normal_vector_one[0]*normal_vector_one[1] - 2*vector_b[0]*vector_b[2]*normal_vector_one[0]*normal_vector_one[2] - vector_b[1]*vector_b[1]*costhetaabc*costhetaabc + vector_b[1]*vector_b[1]*normal_vector_one[0]*normal_vector_one[0] + vector_b[1]*vector_b[1]*normal_vector_one[2]*normal_vector_one[2] - 2*vector_b[1]*vector_b[2]*normal_vector_one[1]*normal_vector_one[2] - vector_b[2]*vector_b[2]*costhetaabc*costhetaabc + vector_b[2]*vector_b[2]*normal_vector_one[0]*normal_vector_one[0] + vector_b[2]*vector_b[2]*normal_vector_one[1]*normal_vector_one[1];
//these two variables appear to be the same. After you get this first version with both variables working, try eliminating one of them and see if the result still works.


    if(fabs(sqrtterm1) < zerotol)
    {
        sqrtterm1 = 0;
    }

    if(fabs(sqrtterm2) < zerotol)
    {
        sqrtterm2 = 0;
    }

    n22 = neg * (vector_b[0]*normal_vector_one[1]*sqrt(sqrtterm1) - vector_b[1]*normal_vector_one[0] * sqrt(sqrtterm2) + neg*(vector_b[0]*vector_b[0]*costhetaabc*normal_vector_one[2] + vector_b[1]*vector_b[1]*costhetaabc*normal_vector_one[2] - vector_b[0]*vector_b[2]*costhetaabc*normal_vector_one[0] - vector_b[1]*vector_b[2]*costhetaabc*normal_vector_one[1])) / (vector_b[0]*vector_b[0]*normal_vector_one[1]*normal_vector_one[1] + vector_b[0]*vector_b[0]*normal_vector_one[2]*normal_vector_one[2] - 2*vector_b[0]*vector_b[1]*normal_vector_one[0]*normal_vector_one[1] - 2*vector_b[0]*vector_b[2]*normal_vector_one[0]*normal_vector_one[2] + vector_b[1]*vector_b[1]*normal_vector_one[0]*normal_vector_one[0] + vector_b[1]*vector_b[1]*normal_vector_one[2]*normal_vector_one[2] - 2*vector_b[1]*vector_b[2]*normal_vector_one[1]*normal_vector_one[2] + vector_b[2]*vector_b[2]*normal_vector_one[0]*normal_vector_one[0] + vector_b[2]*vector_b[2]*normal_vector_one[1]*normal_vector_one[1]);

    n20 = -(vector_b[1]*costhetaabc - vector_b[1]*normal_vector_one[2]*n22 + vector_b[2]*normal_vector_one[1]*n22) / (vector_b[0]*normal_vector_one[1] - vector_b[1]*normal_vector_one[0]);
    n21 = (vector_b[0]*costhetaabc - vector_b[0]*normal_vector_one[2]*n22 + vector_b[2]*normal_vector_one[0]*n22) / (vector_b[0]*normal_vector_one[1] - vector_b[1]*normal_vector_one[0]);

    normal_vector_two[0] = n20;
    normal_vector_two[1] = n21;
    normal_vector_two[2] = n22;

    double planex[3] = {vector_b[0], vector_b[1], vector_b[2]};
    unitize(planex);

    double* planey = cross(normal_vector_two, planex);
    unitize(planey);

    double thetabc = bond_angle*degtorad;

    double costhetabc = cos(thetabc);
    double sinthetabc = sin(thetabc);

    double ux[3] = {costhetabc*planex[0], costhetabc*planex[1], costhetabc*planex[2]};
    double uy[3] = {sinthetabc*planey[0], sinthetabc*planey[1], sinthetabc*planey[2]};

    double u[3] = {ux[0]+uy[0], ux[1]+uy[1], ux[2]+uy[2]};
    unitize(u);

    double newlength = bond_length;

    double v[3] = {newlength*u[0], newlength*u[1], newlength*u[2]};

    point_four[0] = point_one[0] + v[0];
    point_four[1] = point_one[1] + v[1];
    point_four[2] = point_one[2] + v[2];


    return point_four;
}


double placement_checker_POSCAR(POSCAR_struct* POSCAR, molecule_coords* molec_info, int atom_number, double min_scale, double max_scale, double position_and_info[], int z_min, int z_max)
{
    double success = 0; //==1 for succesful placement, ==0 for failure to place succesfully. Automatically set to 0 for safety.
    int a,b,c,d,e; //counter variables for loops

    int terminate = 0; //stops distance checking loop if a collision is found

    double distance1 = 0; //same z neighbors
    double distance2 = 0;
    double distance3 = 0;
    double distance4 = 0;
    double distance5 = 0;
    double distance6 = 0;
    double distance7 = 0;
    double distance8 = 0;
    double distance9 = 0;

    double distance10 = 0; //upstairs z neighbors
    double distance11 = 0;
    double distance12 = 0;
    double distance13 = 0;
    double distance14 = 0;
    double distance15 = 0;
    double distance16 = 0;
    double distance17 = 0;
    double distance18 = 0;

    double distance19 = 0; //downstairs z neighbors
    double distance20 = 0;
    double distance21 = 0;
    double distance22 = 0;
    double distance23 = 0;
    double distance24 = 0;
    double distance25 = 0;
    double distance26 = 0;
    double distance27 = 0;

    int octahedron[8] = {0};

    int close_enough = 0;

    int num_atoms_per_molec = 0;
    int current_atom_index = 1; //counts the current atom's position in any appropriate arrays. Need to subtract 1 for use as an index.

    int actual_number_placed = molec_info->actual_number_placed; //number of molecules of the current type that have been succesfully placed

    double* POSCAR_pos = frac_to_cartesian(POSCAR->latt_vec, POSCAR->par_pos, POSCAR->par_count,POSCAR->scale_factor); //converting fractional lattice coordinates of the POSCAR to Cartesian coordinates.

    //Calculating sm_min and sm_max from the covalent and vdw radii
    int placed_atom_radius_index = return_radius_index(molec_info->structure->atom_ident[atom_number-1]);
    int surface_atom_radius_index;
    int par_count_holder_previous;
    int par_count_holder_current;

    double sm_min;
    double sm_max;

    double l1 = POSCAR->scale_factor*POSCAR->latt_vec[0]; //these are used to create the periodic neighbors' positions
    double l2 = POSCAR->scale_factor*POSCAR->latt_vec[1];
    double l3 = POSCAR->scale_factor*POSCAR->latt_vec[2];
    double l4 = POSCAR->scale_factor*POSCAR->latt_vec[3];
    double l5 = POSCAR->scale_factor*POSCAR->latt_vec[4];
    double l6 = POSCAR->scale_factor*POSCAR->latt_vec[5];
    double l7 = POSCAR->scale_factor*POSCAR->latt_vec[6];
    double l8 = POSCAR->scale_factor*POSCAR->latt_vec[7];
    double l9 = POSCAR->scale_factor*POSCAR->latt_vec[8];

    int holder = 0;
/*
    FILE* cartesian_surface = NULL;

    cartesian_surface = fopen("D://cartesian_surface.txt", "w");

    if(cartesian_surface==NULL)
    {
        printf("\nFailed to open cartesian_surface.txt!");
        exit(123456789);
    }


    for(b=0; b<POSCAR->par_count; b++)
    {
        fprintf(cartesian_surface, "%0.16lf\t%0.16lf\t%0.16lf\n", POSCAR_pos[3*b], POSCAR_pos[3*b+1], POSCAR_pos[3*b+2]);
    }
fclose(cartesian_surface);
*/
/*
    FILE* records = NULL;

    records = fopen("D://records_sm.txt","a");

    if(records==NULL)
    {
        printf("\nFailed to open records_sm.txt!");
        exit(1933);
    }
*/

     molec_info->close_enough = 0;


    for(a=0; a<molec_info->num_atom_types; a++) //Counts the number of atoms that make up the current type of molecule
    {
        num_atoms_per_molec += molec_info->num_each_type[a];
    }


    current_atom_index = (molec_info->actual_number_placed)*(num_atoms_per_molec)+atom_number; //gives the overall atom number of the current atom in the current instance of the current molecule type.
                                                                                      // Ex: H20. If the order is O,H,H, then the first hydrogen on the second molecule would have a current_atom_index of 5 (O,H,H,O,H)
    b=1;
    while((b<=POSCAR->par_count) && (terminate==0))
    { c=0;

        par_count_holder_current = 0;
        par_count_holder_previous = 0;

        for(c=0; c<POSCAR->num_types; c++)
        {
            par_count_holder_current+=POSCAR->par_num[c];

            if((b<=par_count_holder_current) && (b>par_count_holder_previous))
            {
                //printf("\nPOSCAR->par_ident[%i] = %s\tb = %i", c, POSCAR->par_ident[c], b);

                surface_atom_radius_index = return_radius_index(POSCAR->par_ident[c]);
                holder = c;
            }

            par_count_holder_previous = par_count_holder_current;
        }

        sm_min = min_scale*(covalent_radii[placed_atom_radius_index] + covalent_radii[surface_atom_radius_index]);
        sm_max = max_scale*(vdw_radii[placed_atom_radius_index] + vdw_radii[surface_atom_radius_index]);
        //printf("\nsm_min = %lf\tsm_max = %lf", sm_min, sm_max);
        //printf("\nmin_scale = %lf\tmax_scale = %lf", min_scale, max_scale);
        //printf("\nplaced_atom = %s\tsurface_atom = %s", element_ID[placed_atom_radius_index], element_ID[surface_atom_radius_index]);
        //printf("\nplaced_atom covalent = %lf\tsurface_atom covalent = %lf", covalent_radii[placed_atom_radius_index], covalent_radii[surface_atom_radius_index]);
        //printf("\nplaced_atom vdW = %lf\tsurface_atom vdW = %lf\n", vdw_radii[placed_atom_radius_index], vdw_radii[surface_atom_radius_index]);

        //printf("\nPOSCAR->par_ident = %s", POSCAR->par_ident[holder]);

        distance1 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1+l4), (POSCAR_pos[(3*b)-2]-l2+l5), (POSCAR_pos[(3*b)-1]-l3+l6)); //-x, +y

        distance2 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l4), (POSCAR_pos[(3*b)-2]+l5), (POSCAR_pos[(3*b)-1]+l6)); //+y

        distance3 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1+l4), (POSCAR_pos[(3*b)-2]+l2+l5), (POSCAR_pos[(3*b)-1]+l3+l6)); //+x, +y

        distance4 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1), (POSCAR_pos[(3*b)-2]-l2), (POSCAR_pos[(3*b)-1]-l3)); //-x

        distance5 = dist(position_and_info[1], position_and_info[2], position_and_info[3], POSCAR_pos[(3*b)-3], POSCAR_pos[(3*b)-2], POSCAR_pos[(3*b)-1]); //Original unit cell

        distance6 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1), (POSCAR_pos[(3*b)-2]+l2), (POSCAR_pos[(3*b)-1]+l3)); //+x

        distance7 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1-l4), (POSCAR_pos[(3*b)-2]-l2-l5), (POSCAR_pos[(3*b)-1]-l3-l6)); //-x, -y

        distance8 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l4), (POSCAR_pos[(3*b)-2]-l5), (POSCAR_pos[(3*b)-1]-l6)); //-y

        distance9 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1-l4), (POSCAR_pos[(3*b)-2]+l2-l5), (POSCAR_pos[(3*b)-1]+l3-l6)); //+x, -y

    //upstairs z neighbors
        distance10 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1+l4+l7), (POSCAR_pos[(3*b)-2]-l2+l5+l8), (POSCAR_pos[(3*b)-1]-l3+l6+l9)); //-x, +y, +z

        distance11 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l4+l7), (POSCAR_pos[(3*b)-2]+l5+l8), (POSCAR_pos[(3*b)-1]+l6+l9)); //+y, +z

        distance12 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1+l4+l7), (POSCAR_pos[(3*b)-2]+l2+l5+l8), (POSCAR_pos[(3*b)-1]+l3+l6+l9)); //+x, +y, +z

        distance13 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1+l7), (POSCAR_pos[(3*b)-2]-l2+l8), (POSCAR_pos[(3*b)-1]-l3+l9)); //-x, +z

        distance14 = dist(position_and_info[1], position_and_info[2], position_and_info[3], POSCAR_pos[(3*b)-3]+l7, POSCAR_pos[(3*b)-2]+l8, POSCAR_pos[(3*b)-1]+l9); //Original +z

        distance15 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1+l7), (POSCAR_pos[(3*b)-2]+l2+l8), (POSCAR_pos[(3*b)-1]+l3+l9)); //+x, +z

        distance16 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1-l4+l7), (POSCAR_pos[(3*b)-2]-l2-l5+l8), (POSCAR_pos[(3*b)-1]-l3-l6+l9)); //-x, -y, +z

        distance17 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l4+l7), (POSCAR_pos[(3*b)-2]-l5+l8), (POSCAR_pos[(3*b)-1]-l6+l9)); //-y, +z

        distance18 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1-l4+l7), (POSCAR_pos[(3*b)-2]+l2-l5+l8), (POSCAR_pos[(3*b)-1]+l3-l6+l9)); //+x, -y, +z

    //downstairs z neighbors
        distance19 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1+l4-l7), (POSCAR_pos[(3*b)-2]-l2+l5-l8), (POSCAR_pos[(3*b)-1]-l3+l6-l9)); //-x, +y, -z

        distance20 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l4-l7), (POSCAR_pos[(3*b)-2]+l5-l8), (POSCAR_pos[(3*b)-1]+l6-l9)); //+y, -z

        distance21 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1+l4-l7), (POSCAR_pos[(3*b)-2]+l2+l5-l8), (POSCAR_pos[(3*b)-1]+l3+l6-l9)); //+x, +y, -z

        distance22 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1-l7), (POSCAR_pos[(3*b)-2]-l2-l8), (POSCAR_pos[(3*b)-1]-l3-l9)); //-x, -z

        distance23 = dist(position_and_info[1], position_and_info[2], position_and_info[3], POSCAR_pos[(3*b)-3]-l7, POSCAR_pos[(3*b)-2]-l8, POSCAR_pos[(3*b)-1]-l9); //Original -z

        distance24 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1-l7), (POSCAR_pos[(3*b)-2]+l2-l8), (POSCAR_pos[(3*b)-1]+l3-l9)); //+x, -z

        distance25 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l1-l4-l7), (POSCAR_pos[(3*b)-2]-l2-l5-l8), (POSCAR_pos[(3*b)-1]-l3-l6-l9)); //-x, -y, -z

        distance26 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]-l4-l7), (POSCAR_pos[(3*b)-2]-l5-l8), (POSCAR_pos[(3*b)-1]-l6-l9)); //-y, -z

        distance27 = dist(position_and_info[1], position_and_info[2], position_and_info[3], (POSCAR_pos[(3*b)-3]+l1-l4-l7), (POSCAR_pos[(3*b)-2]+l2-l5-l8), (POSCAR_pos[(3*b)-1]+l3-l6-l9)); //+x, -y, -z


        if(position_and_info[3]<z_min)
        {
            terminate = 1;
        }

        if(position_and_info[3]>z_max)
        {
            terminate = 1;
        }


        if(distance5<sm_min) //this if-statement is separate from the following one so that different error messages may be easily printed for the different scenarios
        {
            terminate = 1;
//printf("\nfailure 1, distance5 = %lf, sm_min = %lf, atom_number = %i", distance5, sm_min, atom_number);

//printf("\nPOSCAR_x = %lf\tPOSCAR_y = %lf\tPOSCAR_z = %lf", POSCAR_pos[3*b-3], POSCAR_pos[3*b-2], POSCAR_pos[3*b-1]);
            if(molec_info->replacing==1)
            {
                for(e=0; e<(molec_info->structure->num_atoms); e++)
                { //printf("\nb = %i\treplaced_par_pos_indices[%i] = %i", b, e, molec_info->replaced_par_pos_indices[e]);
                    if((b-1)==molec_info->replaced_par_pos_indices[e])
                    { //printf("\nA replaced particle is overlapping with its original copy1");

                        terminate = 0;
                    }
                }
            }
            //printf("\tToo close to surface");
        }
            //checking against neighboring periodic copies of surface
        if(((distance1<sm_min) || (distance2<sm_min) || (distance3<sm_min) || (distance4<sm_min) || (distance6<sm_min) || (distance7<sm_min) || (distance8<sm_min) || (distance9<sm_min)) && (terminate==0))
        {
            terminate = 1;
//printf("\nfailure 2");
            if(molec_info->replacing==1)
            {
                for(e=0; e<(molec_info->structure->num_atoms); e++)
                {
                    if((b-1)==molec_info->replaced_par_pos_indices[e])
                    { //printf("\tA replaced particle is overlapping with its original copy2");
                        terminate = 0;
                    }
                }
            }
            //printf("\tToo close to periodic copy/copies of original surface unit cell");
        }

        if(((distance10<sm_min) || (distance11<sm_min) || (distance12<sm_min) || (distance13<sm_min) || (distance14<sm_min) || (distance15<sm_min) || (distance16<sm_min) || (distance17<sm_min) || (distance18<sm_min)) && (terminate==0))
        {
            terminate = 1;
//printf("\nfailure 3");
            if(molec_info->replacing==1)
            {
                for(e=0; e<(molec_info->structure->num_atoms); e++)
                {
                    if((b-1)==molec_info->replaced_par_pos_indices[e])
                    { //printf("\tA replaced particle is overlapping with its original copy3");
                        terminate = 0;
                    }
                }
            }
            //printf("\tToo close to periodic copy/copies of original surface unit cell");
        }

        if(((distance19<sm_min) || (distance20<sm_min) || (distance21<sm_min) || (distance22<sm_min) || (distance23<sm_min) || (distance24<sm_min) || (distance25<sm_min) || (distance26<sm_min) || (distance27<sm_min)) && (terminate==0))
        {
            terminate = 1;
//printf("\nfailure 4");
            if(molec_info->replacing==1)
            {
                for(e=0; e<(molec_info->structure->num_atoms); e++)
                {
                    if((b-1)==molec_info->replaced_par_pos_indices[e])
                    { //printf("\tA replaced particle is overlapping with its original copy4");
                        terminate = 0;
                    }
                }
            }
            //printf("\tToo close to periodic copy/copies of original surface unit cell");
        }

        if((distance1<sm_max) || (distance2<sm_max) || (distance3<sm_max) || (distance4<sm_max) || (distance5<sm_max) || (distance6<sm_max) || (distance7<sm_max) || (distance8<sm_max) || (distance9<sm_max))
        {
            molec_info->close_enough = 1;
        }
/*
        if(molec_info->replacing==1)
        {
            fprintf(records, "\n%0.16lf\t%0.16lf\t%0.16lf", position_and_info[1], position_and_info[2], position_and_info[3]); //New x, new y, newz
            fprintf(records, "\t%0.16lf\t%0.16lf\t%0.16lf", POSCAR_pos[(3*b)-3], POSCAR_pos[(3*b)-2], POSCAR_pos[(3*b)-1]); // POSCAR x, POSCAR y, POSCAR z
            fprintf(records, "\t%lf\t%lf", sm_min, sm_max);
            fprintf(records, "\t%lf", distance5);
            fprintf(records, "\t%s\t%s", element_ID[placed_atom_radius_index], element_ID[surface_atom_radius_index]);
            fprintf(records, "\t%i", terminate);
        }
*/


        b++;
    } ///end of loop
    if(terminate==0)
    {
        success = 1;
    }

    //printf("\n\nPOSCAR close_enough = %i", molec_info->close_enough);
    //printf("\nPOSCAR checking success = %.0f", success);

    //fclose(records);
    free(POSCAR_pos);

    return success;
}


double placement_checker_placedparticles(POSCAR_struct* POSCAR, placed_particles* par_list, molecule_coords* molec_info, int atom_number, double min_scale, double max_scale, double position_and_info[])
{
/*
    FILE* fileptr_oxygen = NULL;

    fileptr_oxygen = fopen("D://oxygen_real_time.txt", "a");

    if(fileptr_oxygen==NULL)
    {
        printf("\nFailed to open oxygen_real_time.txt!");
        exit(1738);
    }

    FILE* fileptr_hydrogen = NULL;

    fileptr_hydrogen = fopen("D://hydrogen_real_time.txt", "a");

    if(fileptr_hydrogen==NULL)
    {
        printf("\nFailed to open hydrogen_real_time.txt!");
        exit(1738);
    }
*/



//////////////////////////////////////////

    double success = 1;

    double mm_min = 1;
    double mm_max = 1;

    int a,b,c; //counter variables for loops

    int num_types_molecules = par_list->num_types_molecules;
    int num_atoms_per_molecule = 0;
    int num_atoms_placed = 0; ///IMPORTANT- This is the number of atoms in the molecules that have been succesfully and completely placed,
                                    /// and does NOT include the atoms in the molecule currently being placed.

    double distance1 = 0; //same z neighbors
    double distance2 = 0;
    double distance3 = 0;
    double distance4 = 0;
    double distance5 = 0;
    double distance6 = 0;
    double distance7 = 0;
    double distance8 = 0;
    double distance9 = 0;

    double distance10 = 0; //upstairs z neighbors
    double distance11 = 0;
    double distance12 = 0;
    double distance13 = 0;
    double distance14 = 0;
    double distance15 = 0;
    double distance16 = 0;
    double distance17 = 0;
    double distance18 = 0;

    double distance19 = 0; //downstairs z neighbors
    double distance20 = 0;
    double distance21 = 0;
    double distance22 = 0;
    double distance23 = 0;
    double distance24 = 0;
    double distance25 = 0;
    double distance26 = 0;
    double distance27 = 0;

    double l1 = POSCAR->scale_factor*POSCAR->latt_vec[0]; //these are used to create the periodic neighbors' positions
    double l2 = POSCAR->scale_factor*POSCAR->latt_vec[1];
    double l3 = POSCAR->scale_factor*POSCAR->latt_vec[2];
    double l4 = POSCAR->scale_factor*POSCAR->latt_vec[3];
    double l5 = POSCAR->scale_factor*POSCAR->latt_vec[4];
    double l6 = POSCAR->scale_factor*POSCAR->latt_vec[5];
    double l7 = POSCAR->scale_factor*POSCAR->latt_vec[6];
    double l8 = POSCAR->scale_factor*POSCAR->latt_vec[7];
    double l9 = POSCAR->scale_factor*POSCAR->latt_vec[8];

    int num_comparisons = 0;

    int override = 0;
    int terminate = 0;
    int too_close = 0;

    int mobile_particle_radius_index = 0; //compared_to
    int placed_particle_radius_index = 0; //currently being placed

    double box_x = POSCAR->scale_factor*(POSCAR->latt_vec[0] + POSCAR->latt_vec[3] + POSCAR->latt_vec[6]);
    double box_y = POSCAR->scale_factor*(POSCAR->latt_vec[1] + POSCAR->latt_vec[4] + POSCAR->latt_vec[7]);


    a=0;
    while((a<num_types_molecules) && (terminate==0))
    {

    num_atoms_per_molecule = 0;
        for(b=0; b<par_list->molecule_info[a].num_atom_types; b++)
        {
            num_atoms_per_molecule += par_list->molecule_info[a].num_each_type[b]; //calculates the number of atoms in a molecule of the current type
        }

        num_atoms_placed = (num_atoms_per_molecule * par_list->molecule_info[a].actual_number_placed); //the total number of atoms belonging to the current molecular type that have been placed

        if(num_atoms_placed==0)
        {
            terminate = 1;
            override = 1;
        }
        c = 1;
    num_comparisons = 0;
        while((c<=num_atoms_placed) && (terminate==0))
        { num_comparisons++;

            mobile_particle_radius_index = c;

            while(mobile_particle_radius_index>par_list->molecule_info[a].structure->num_atoms)
            {
                mobile_particle_radius_index-=par_list->molecule_info[a].structure->num_atoms;
            }
//printf("\n\nmpri = %i\ta = %i\tc= %i\tnum_atoms = %i", mobile_particle_radius_index, a, c, par_list->molecule_info[a].structure->num_atoms);
            //printf("\nplacing %s", molec_info->structure->atom_ident[atom_number-1]);
            //printf("\ncomparing to %s", par_list->molecule_info[a].structure->atom_ident[mobile_particle_radius_index-1]);


            mobile_particle_radius_index = return_radius_index(par_list->molecule_info[a].structure->atom_ident[mobile_particle_radius_index-1]); //printf("\n2577- mpri = %i", mobile_particle_radius_index);
            placed_particle_radius_index = return_radius_index(molec_info->structure->atom_ident[atom_number-1]); //printf("\n2578\t ppri = %i\tmpri = %i", placed_particle_radius_index, mobile_particle_radius_index);
            //printf("\nmobile_particle_radius_index = %i", mobile_particle_radius_index);
            //printf("\ncompared_to ident = %s", molec_info->structure->atom_ident[atom_number-1]);
            mm_min = min_scale*(covalent_radii[placed_particle_radius_index] + covalent_radii[mobile_particle_radius_index]);
            mm_max = max_scale*(vdw_radii[placed_particle_radius_index] + vdw_radii[mobile_particle_radius_index]);


            distance1 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1+l4, par_list->molecule_info[a].coords[(3*c)-2]-l2+l5, par_list->molecule_info[a].coords[(3*c)-1]-l3+l6); //-x, +y

            distance2 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l4, par_list->molecule_info[a].coords[(3*c)-2]+l5, par_list->molecule_info[a].coords[(3*c)-1]+l6); //+y

            distance3 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1+l4, par_list->molecule_info[a].coords[(3*c)-2]+l2+l5, par_list->molecule_info[a].coords[(3*c)-1]+l3+l6); //+x, +y

            distance4 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1, par_list->molecule_info[a].coords[(3*c)-2]-l2, par_list->molecule_info[a].coords[(3*c)-1]-l3); //-x

            distance5 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3], par_list->molecule_info[a].coords[(3*c)-2], par_list->molecule_info[a].coords[(3*c)-1]); //Original unit cell

            distance6 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1, par_list->molecule_info[a].coords[(3*c)-2]+l2, par_list->molecule_info[a].coords[(3*c)-1]+l3); //+x

            distance7 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1-l4, par_list->molecule_info[a].coords[(3*c)-2]-l2-l5, par_list->molecule_info[a].coords[(3*c)-1]-l3-l6); //-x, -y

            distance8 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l4, par_list->molecule_info[a].coords[(3*c)-2]-l5, par_list->molecule_info[a].coords[(3*c)-1]-l6); //-y

            distance9 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1-l4, par_list->molecule_info[a].coords[(3*c)-2]+l2-l5, par_list->molecule_info[a].coords[(3*c)-1]+l3-l6); //+x, -y


        //upstairs z neighbors
            distance10 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1+l4+l7, par_list->molecule_info[a].coords[(3*c)-2]-l2+l5+l8, par_list->molecule_info[a].coords[(3*c)-1]-l3+l6+l9); //-x, +y, +z

            distance11 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l4+l7, par_list->molecule_info[a].coords[(3*c)-2]+l5+l8, par_list->molecule_info[a].coords[(3*c)-1]+l6+l9); //+y, +z

            distance12 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1+l4+l7, par_list->molecule_info[a].coords[(3*c)-2]+l2+l5+l8, par_list->molecule_info[a].coords[(3*c)-1]+l3+l6+l9); //+x, +y, +z

            distance13 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1+l7, par_list->molecule_info[a].coords[(3*c)-2]-l2+l8, par_list->molecule_info[a].coords[(3*c)-1]-l3+l9); //-x, +z

            distance14 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l7, par_list->molecule_info[a].coords[(3*c)-2]+l8, par_list->molecule_info[a].coords[(3*c)-1]+l9); //Original +z

            distance15 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1+l7, par_list->molecule_info[a].coords[(3*c)-2]+l2+l8, par_list->molecule_info[a].coords[(3*c)-1]+l3+l9); //+x, +z

            distance16 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1-l4+l7, par_list->molecule_info[a].coords[(3*c)-2]-l2-l5+l8, par_list->molecule_info[a].coords[(3*c)-1]-l3-l6+l9); //-x, -y, +z

            distance17 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l4+l7, par_list->molecule_info[a].coords[(3*c)-2]-l5+l8, par_list->molecule_info[a].coords[(3*c)-1]-l6+l9); //-y, +z

            distance18 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1-l4+l7, par_list->molecule_info[a].coords[(3*c)-2]+l2-l5+l8, par_list->molecule_info[a].coords[(3*c)-1]+l3-l6+l9); //+x, -y, +z


        //downstairs z neighbors
            distance19 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1+l4-l7, par_list->molecule_info[a].coords[(3*c)-2]-l2+l5-l8, par_list->molecule_info[a].coords[(3*c)-1]-l3+l6-l9); //-x, +y, -z

            distance20 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l4-l7, par_list->molecule_info[a].coords[(3*c)-2]+l5-l8, par_list->molecule_info[a].coords[(3*c)-1]+l6-l9); //+y, -z

            distance21 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1+l4-l7, par_list->molecule_info[a].coords[(3*c)-2]+l2+l5-l8, par_list->molecule_info[a].coords[(3*c)-1]+l3+l6-l9); //+x, +y, -z

            distance22 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1-l7, par_list->molecule_info[a].coords[(3*c)-2]-l2-l8, par_list->molecule_info[a].coords[(3*c)-1]-l3-l9); //-x, -z

            distance23 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l7, par_list->molecule_info[a].coords[(3*c)-2]-l8, par_list->molecule_info[a].coords[(3*c)-1]-l9); //Original -z

            distance24 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1-l7, par_list->molecule_info[a].coords[(3*c)-2]+l2-l8, par_list->molecule_info[a].coords[(3*c)-1]+l3-l9); //+x, -z

            distance25 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l1-l4-l7, par_list->molecule_info[a].coords[(3*c)-2]-l2-l5-l8, par_list->molecule_info[a].coords[(3*c)-1]-l3-l6-l9); //-x, -y, -z

            distance26 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]-l4-l7, par_list->molecule_info[a].coords[(3*c)-2]-l5-l8, par_list->molecule_info[a].coords[(3*c)-1]-l6-l9); //-y, -z

            distance27 = dist(position_and_info[1], position_and_info[2], position_and_info[3], par_list->molecule_info[a].coords[(3*c)-3]+l1-l4-l7, par_list->molecule_info[a].coords[(3*c)-2]+l2-l5-l8, par_list->molecule_info[a].coords[(3*c)-1]+l3-l6-l9); //+x, -y, -z


            if((distance1<mm_min) || (distance2<mm_min) || (distance3<mm_min) || (distance4<mm_min) || (distance5<mm_min) || (distance6<mm_min) || (distance7<mm_min) || (distance8<mm_min) || (distance9<mm_min))
            { //printf("\n\n\t\tmp-mp distance5 = %lf\n\n", distance5);
                //printf("\tMP overlap\n");
                too_close = 1;
                terminate = 1;
                success = 0;
                override = 0;
            }

            if((distance10<mm_min) || (distance11<mm_min) || (distance12<mm_min) || (distance13<mm_min) || (distance14<mm_min) || (distance15<mm_min) || (distance16<mm_min) || (distance17<mm_min) || (distance18<mm_min))
            { //printf("\n\n\t\tmp-mp distance5 = %lf\n\n", distance5);
                //printf("\tMP overlap\n");
                too_close = 1;
                terminate = 1;
                success = 0;
                override = 0;
            }

            if((distance19<mm_min) || (distance20<mm_min) || (distance21<mm_min) || (distance22<mm_min) || (distance23<mm_min) || (distance24<mm_min) || (distance25<mm_min) || (distance26<mm_min) || (distance27<mm_min))
            { //printf("\n\n\t\tmp-mp distance5 = %lf\n\n", distance5);
                //printf("\tMP overlap\n");
                too_close = 1;
                terminate = 1;
                success = 0;
                override = 0;
            }

            if(molec_info->close_enough==0)
            {
                if((distance1<mm_max) || (distance2<mm_max) || (distance3<mm_max) || (distance4<mm_max) || (distance5<mm_max) || (distance6<mm_max) || (distance7<mm_max) || (distance8<mm_max) || (distance9<mm_max))
                {
                    molec_info->close_enough = 1;
                }
            }
            c++;
        } ///end of loop

    a++;
    } ///end of loop



    if((terminate==1) || (molec_info->close_enough==0))
    {
        success = 0;
    }

    if((terminate==1) && (override==1) && (molec_info->close_enough==1))
    {
        success = 1;
    }

    if(too_close==1)
    {
        success = 0;
    }
/*
    if(success==1)
    {
        printf("\nsuccess = %.0f\t\tclose_enough = %i", success, molec_info->close_enough);
    }


    if(success==1)
    {
        if(atom_number==1)
        {
            fprintf(fileptr_oxygen, "%0.16lf   %0.16lf   %0.16lf\n", position_and_info[1], position_and_info[2], position_and_info[3]);
        }
        else if((atom_number==2) || (atom_number==3))
        {
            fprintf(fileptr_hydrogen, "%0.16lf   %0.16lf   %0.16lf\n", position_and_info[1], position_and_info[2], position_and_info[3]);
        }
    }
*/

    //printf("\nPlaced particles checking success = %.0f", success);
    //fclose(records);
    //fclose(fileptr_hydrogen);
    //fclose(fileptr_oxygen);
    return success;
}


void fill_actual_number_placed(placed_particles* par_list, POSCAR_struct* resubmission)
{
    int a, b;
    int actual_number_placed;
    int num_atoms;

    int molecule_count = 0;

    for(a=0; a<par_list->num_types_molecules; a++)
    {
        actual_number_placed = 0;
        num_atoms = par_list->molecule_info[a].structure->num_atoms;

        for(b=0; b<resubmission->par_count; b++)
        {
            if(resubmission->molecular_type_tag[b]==a)
            {
                actual_number_placed++;
            }
        }
        par_list->molecule_info[a].actual_number_placed = actual_number_placed/num_atoms;

        molecule_count+= par_list->molecule_info[a].actual_number_placed;
    }

    par_list->molecule_count = molecule_count;

}

  ///////////////////////////////////////////////
 ////    Output and organizing functions    ////
///////////////////////////////////////////////

int organize_atom_idents(placed_particles* particle_list, POSCAR_struct* input)
{

    /*

        The purpose of this function is to consolidate the POSCAR and placed particles information in order to determine:
            a) How many different elements are in the system
            b) The user-specified maximum number of each element (not directly specified, but calculated from maximum molecule numbers given)
            c) Begin organizing information for POSCAR writing, which will be done after all molecules have been placed about the original POSCAR.

    */
    int num_types_molecules = particle_list->num_types_molecules;
    int num_types_atoms = 0;

    int num_atom_types; //holds the number of types of different atoms in each molecule. The original value is contained in the molecule_coords struct
    int poscar_num_types = input->num_types;
    int a;

    int current_molec_type; //cycles through all molecule types
    int current_atom_type; //cycles through atoms in a particular molecule type

    int current_ident;


    int cmp_result = 0;


/////////////////////////////////////////////////////////////////////////////

    //Adds POSCAR information to num_types_atoms, master_max_num_each_type, and master_ident in placed_particles structure

    for(a=0; a<input->num_types; a++)
    {
        particle_list->master_max_num_each_type[a] = input->par_num[a];

        copyString(particle_list->master_ident[a], input->par_ident[a]);
    }

    particle_list->master_num_types_atoms = input->num_types;

    for(a=input->num_types; a<112; a++)
    {
        particle_list->master_max_num_each_type[a] = 0;
    }

////////////////////////////////////////////////////////////////////////////

    //Adds molecule information to the aforementioned data things

    for(current_molec_type=0; current_molec_type<num_types_molecules; current_molec_type++)
    {


        num_types_atoms = particle_list->molecule_info[current_molec_type].num_atom_types;

        for(current_atom_type=0; current_atom_type<(particle_list->molecule_info[current_molec_type].num_atom_types); current_atom_type++)
        {   //printf("\nnum_atom_types = %i", particle_list->molecule_info[current_molec_type]->num_atom_types);
            //printf("\ncurrent atom type = %s\t%i", particle_list->molecule_info[current_molec_type].type_ident[current_atom_type], current_atom_type);

            for(current_ident=0; current_ident<particle_list->master_num_types_atoms; current_ident++)
            {
                cmp_result = strcmp(particle_list->master_ident[current_ident], particle_list->molecule_info[current_molec_type].type_ident[current_atom_type]);


                if(cmp_result==0)
                { //printf("\n\t\t\tMatch");
                    particle_list->master_max_num_each_type[current_ident]+=particle_list->molecule_info[current_molec_type].max_num_each_type[current_atom_type];
                    current_ident = particle_list->master_num_types_atoms;
                }
                else if(cmp_result!=0)
                { //printf("\n\t\t\tMismatch");

                    if(current_ident==(particle_list->master_num_types_atoms-1))
                    {
                        //printf("\n\tCreating new type in master array");
                        //printf("\n\tcurrent_ident = %i\tmaster_num_types_atoms = %i", current_ident, particle_list->master_num_types_atoms);
                        copyString(particle_list->master_ident[particle_list->master_num_types_atoms], particle_list->molecule_info[current_molec_type].type_ident[current_atom_type]);
                        particle_list->master_max_num_each_type[particle_list->master_num_types_atoms] += particle_list->molecule_info[current_molec_type].max_num_each_type[current_atom_type];
                        particle_list->master_num_types_atoms+=1;
                        current_ident++;
                    }
                }
            }
        }
    }
/*
    //printf("\n\n");
    for(a=0; a<particle_list->master_num_types_atoms; a++)
    {
        printf("\n%i\t%s", particle_list->master_max_num_each_type[a], particle_list->master_ident[a]);
    }
*/

    return 0;
}


void POSCAR_element_sorter(POSCAR_struct* POSCAR, placed_particles* par_list)
{

    int a, b, c, d; //counter variables

    int terminate_inner = 0; //ends inner sorting loop if finished before loop condition(s) satisfied
    int terminate_middle = 0; //ends middle sorting loop if finished before loop condition(s) satisfied
    //Go back through these loops and incorporate terminate_xxxx functionalities into loop conditions.

    par_list->master_actual_num_each_type = (int*) malloc(sizeof(int)*par_list->master_num_types_atoms);

    for(a=0; a<par_list->master_num_types_atoms; a++)
    {
        par_list->master_actual_num_each_type[a] = 0;
        //printf("\n\t%s", par_list->master_ident[a]);
    }

    for(a=0; a<POSCAR->num_types; a++)
    {
        par_list->master_actual_num_each_type[a] = POSCAR->par_num[a];
        //printf("\n\n\nelement = %s", par_list->master_ident[a]);
        //printf("\n\tnum = %i\n", par_list->master_actual_num_each_type[a]);
    }



c=0;
d=0;

//printf("\npar_list->num_types_molecules = %i", par_list->num_types_molecules);
//printf("\npar_list->master_num_types_atoms = %i", par_list->master_num_types_atoms);
//printf("\npar_list->molecule_info[0]->num_atom_types = %i", par_list->molecule_info[0]->num_atom_types);
    for(b=0; b<par_list->num_types_molecules; b++)
    {
        while((c<par_list->master_num_types_atoms) && (terminate_middle==0))
        {
            while((d<par_list->molecule_info[b].num_atom_types) && (terminate_inner==0))
            {
                //printf("\n\ttype 1 = %s\ttype 2 = %s", par_list->master_ident[c], par_list->molecule_info[b]->type_ident[d]);

                if(strcmp(par_list->master_ident[c], par_list->molecule_info[b].type_ident[d])==0)
                {
                    //printf("\n\t\tMatch\n");
                    //printf("\nAdded %i to master %s atom count", (par_list->molecule_info[b].num_each_type[d])*(par_list->molecule_info[b].actual_number_placed),  par_list->master_ident[c]);
                    //printf("\nactual_number_placed = %i", par_list->molecule_info[b].actual_number_placed);
                    //printf("\nnum %s atoms per molec = %i", par_list->master_ident[c], par_list->molecule_info[b].num_each_type[d]);
                    par_list->master_actual_num_each_type[c] += (par_list->molecule_info[b].num_each_type[d])*(par_list->molecule_info[b].actual_number_placed);
                    terminate_inner = 1;
                }
                d++;
                terminate_inner = 0;
            }

            d=0;
            c++;
            terminate_middle = 0;
        }

        c=0;
    }
/*
printf("\n\n");

int e;

for(e=0; e<par_list->master_num_types_atoms; e++)
{
    printf("\nmaster_actual_num_each_type[%i] = %i", e, par_list->master_actual_num_each_type[e]);
}
*/

}


void output_sorter(POSCAR_struct* input, placed_particles* par_list, POSCAR_struct* output)
{
    /**This function sorts and moves the coordinates from the coords pointer in par_list->molecule_info over to the output POSCAR**/



    //printf("\n\n________BEGINNING OUTPUT SORTING ROUTINE________\n\n");
    int total_num_atoms = 0;

    int a, b, c, d, current_atom, e, f, g, h, i, j, k, l, m, n, o, p, q ,r, s, t, u;

    total_num_atoms+= input->par_count;

    for(a=0; a<par_list->num_types_molecules; a++)
    {
        total_num_atoms += (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms);
    }

    //printf("\ntotal_num_atoms = %i\n", total_num_atoms);

    int* molecular_ID = (int*) malloc(total_num_atoms*sizeof(int));
    int* molecular_type_tag = (int*) malloc(total_num_atoms*sizeof(int));
    int* atom_ID_tag = (int*) malloc(total_num_atoms*sizeof(int));
    int* bonded_to = (int*) malloc(2*total_num_atoms*sizeof(int));
    int* angle = (int*) malloc(3*total_num_atoms*sizeof(int));

    double* positions = (double*) malloc(3*total_num_atoms*sizeof(double)); //Placed particles' positions are initially Cartesian. Convert them to the original lattice coordinate system
    double* pp_pos = (double*) malloc(3*(total_num_atoms-(input->par_count))*sizeof(double));

    int* atom_ID = (int*) malloc(total_num_atoms*sizeof(int));
    int* POSCAR_ID = (int*) malloc(input->par_count*sizeof(int));
    int* pp_ID = (int*) malloc((total_num_atoms-(input->par_count))*sizeof(int));



    //Filling pp_ID array
    int pp_ID_index = 0;
    int POSCAR_ID_index = 0;

    for(h=0; h<par_list->num_types_molecules; h++) //filling pp_ID
    {
        for(i=0; i<par_list->molecule_info[h].actual_number_placed; i++)
        {
            for(j=0; j<par_list->molecule_info[h].structure->num_atoms; j++)
            {
                for(k=0; k<par_list->master_num_types_atoms; k++)
                { //printf("\nmaster_ident[%i] = %s", k, par_list->master_ident[k]);
                    if(strcmp(par_list->molecule_info[h].structure->atom_ident[j], par_list->master_ident[k])==0)
                    {
                        pp_ID[pp_ID_index] = k;
                        //printf("\nMATCH");
                        //printf("\n\tmaster_ident[%i] = %s\tpp_ID[%i] = %i",k, par_list->master_ident[k], pp_ID_index, pp_ID[pp_ID_index]);
                        pp_ID_index++;
                    }
                }
            }
        }
    }
/*
for(k=0; k<par_list->master_num_types_atoms; k++)
{
    printf("\nmaster_ident[%i] = %s", k, par_list->master_ident[k]);
}

for(l=0; l<pp_ID_index; l++)
{
    printf("\npp_ID[%i] = %i", l, pp_ID[l]);
}
*/

    int* poscar_elements = (int*) malloc(input->par_count*sizeof(int));
    int current_index = 0;

    for(s=0; s<input->par_count; s++)
    {
        if(s<input->par_num[current_index])
        {
            poscar_elements[s] = current_index;
        }
        else
        {
            current_index++;
        }
    }


    for(m=0; m<input->num_types; m++) //filling POSCAR_ID
    {
        for(n=0; n<input->par_num[m]; n++)
        {
            POSCAR_ID[POSCAR_ID_index] = m;
            //printf("\nPOSCAR_ID[%i] = %i", POSCAR_ID_index, POSCAR_ID[POSCAR_ID_index]);
            POSCAR_ID_index++;
        }
    }

/*
    for(l=0; l<input->num_types; l++)
    {
        printf("\npar_ident[%i] = %s", l, par_list->master_ident[l]);
    }
*/

    for(l=0; l<input->par_count; l++) //Filling atom_ID from POSCAR_ID
    {
        atom_ID[l] = POSCAR_ID[l];
    }

    for(o=input->par_count; o<total_num_atoms; o++) //Finishing off atom_ID with pp_ID
    {
        atom_ID[o] = pp_ID[(o-(input->par_count))];
        //printf("\natom_ID[%i] = %i", o, atom_ID[o]);
    }


///////////////////////////////////////////////////////////


    //Filling molecular_ID and position arrays

    for(b=0; b<input->par_count; b++)
    {
        molecular_ID[b] = input->molecular_ID[b];

        molecular_type_tag[b] = -1;

        positions[3*b] = input->par_pos[3*b];
        positions[3*b+1] = input->par_pos[3*b+1];
        positions[3*b+2] = input->par_pos[3*b+2];
    }

    current_atom = input->par_count;

    for(c=0; c<par_list->num_types_molecules; c++)
    {
        for(d=0; d<((par_list->molecule_info[c].actual_number_placed)*(par_list->molecule_info[c].structure->num_atoms)); d++, current_atom++) //this loop assembles the initial versions of the final molecular_ID, position, and atom_ID arrays
        {
            molecular_ID[current_atom] = par_list->molecule_info[c].molecular_ID[d];

            molecular_type_tag[current_atom] = c;

            positions[3*current_atom] = par_list->molecule_info[c].coords[3*d];
            positions[3*current_atom+1] = par_list->molecule_info[c].coords[3*d+1];
            positions[3*current_atom+2] = par_list->molecule_info[c].coords[3*d+2];

        }
    }

    for(f=0; f<(total_num_atoms-(input->par_count)); f++) //this loop moves the placed particle positions to a separate array for coordinate system conversion
    {
        pp_pos[3*f] = positions[3*(f+input->par_count)];
        pp_pos[3*f+1] = positions[3*(f+input->par_count)+1];
        pp_pos[3*f+2] = positions[3*(f+input->par_count)+2];

        //printf("\npp_pos[%i] = %lf", 3*f, pp_pos[3*f]);
        //printf("\npp_pos[%i] = %lf", 3*f+1, pp_pos[3*f+1]);
        //printf("\npp_pos[%i] = %lf", 3*f+2, pp_pos[3*f+2]);
    }

    double* pp_pos_frac = NULL;

    pp_pos_frac = cartesian_to_frac(input->latt_vec, pp_pos, total_num_atoms-input->par_count, input->scale_factor); //coordinate system conversion (Cartesian to original lattice coordinates)

    for(g=0; g<(total_num_atoms-(input->par_count)); g++) //this loop moves the placed particle positions back to the main array after coordinate system conversion
    {
        positions[3*(g+input->par_count)] = pp_pos_frac[3*g];
        positions[3*(g+input->par_count)+1] = pp_pos_frac[3*g+1];
        positions[3*(g+input->par_count)+2] = pp_pos_frac[3*g+2];

        if(positions[3*(g+input->par_count)]<0)
        {
            positions[3*(g+input->par_count)]+=1;
        }
        else if(positions[3*(g+input->par_count)]>1)
        {
            positions[3*(g+input->par_count)]-=1;
        }

        if(positions[3*(g+input->par_count)+1]<0)
        {
            positions[3*(g+input->par_count)+1]+=1;
        }
        else if(positions[3*(g+input->par_count)+1]>1)
        {
            positions[3*(g+input->par_count)+1]-=1;
        }

        if(positions[3*(g+input->par_count)+2]<0)
        {
            positions[3*(g+input->par_count)+2]+=1;
        }
        else if(positions[3*(g+input->par_count)+2]>1)
        {
            positions[3*(g+input->par_count)+2]-=1;
        }

    }


/*
    for(e=0; e<total_num_atoms; e++)
    {
        //printf("\nmolecular_ID[%i] = %i", e, molecular_ID[e]);
        printf("\npositions[%i] = %lf", 3*e, positions[3*e]);
        printf("\npositions[%i] = %lf", 3*e+1, positions[3*e+1]);
        printf("\npositions[%i] = %lf", 3*e+2, positions[3*e+2]);
    }
*/



 ////////   Sorting Output Arrays   ////////


    int* order_of_final_arrays = (int*) malloc(total_num_atoms*sizeof(int));

    output->molecular_type_tag = (int*) malloc(total_num_atoms*sizeof(int));
    output->molecular_ID = (int*) malloc(total_num_atoms*sizeof(int));
    output->par_pos = (double*) malloc(3*total_num_atoms*sizeof(double));
    output->par_num = (int*) malloc(par_list->master_num_types_atoms*sizeof(int));
    output->atom_ID = (int*) malloc(total_num_atoms*sizeof(int));
    output->bonded_to = (int*) malloc(2*total_num_atoms*sizeof(int));
    output->angle = (int*) malloc(3*total_num_atoms*sizeof(int));

    int aa, bb;
    int sorting_index = 0;

    for(aa=0; aa<par_list->master_num_types_atoms; aa++)
    {
        for(bb=0; bb<total_num_atoms; bb++)
        {
            if(atom_ID[bb]==aa)
            {
                order_of_final_arrays[sorting_index] = bb;
                sorting_index++;
            }
        }
    }



  ////////   Filling Output Structure   ////////

    for(aa=0; aa<input->par_count; aa++)
    {
        bonded_to[2*aa] = -1;
        bonded_to[2*aa+1] = -1;

        angle[3*aa] = -1;
        angle[3*aa+1] = -1;
        angle[3*aa+2] = -1;

        atom_ID_tag[aa] = -1;
    }


    int current_pos = input->par_count;

    for(aa=0; aa<par_list->num_types_molecules; aa++)
    {
        for(bb=0; bb<(par_list->molecule_info[aa].structure->num_atoms*par_list->molecule_info[aa].actual_number_placed); bb++)
        {
            bonded_to[2*current_pos] = par_list->molecule_info[aa].bonded_to[2*bb];
            bonded_to[2*current_pos+1] = par_list->molecule_info[aa].bonded_to[2*bb+1];

            angle[3*current_pos] = par_list->molecule_info[aa].angle[3*bb];
            angle[3*current_pos+1] = par_list->molecule_info[aa].angle[3*bb+1];
            angle[3*current_pos+2] = par_list->molecule_info[aa].angle[3*bb+2];

            atom_ID_tag[current_pos] = par_list->molecule_info[aa].atom_ID[bb];

            current_pos++;
        }
    }


    for(aa=0; aa<total_num_atoms; aa++)
    {
        //printf("\norder_of_final_arrays[%i] = %i", aa, order_of_final_arrays[aa]);

        output->molecular_ID[aa] = molecular_ID[order_of_final_arrays[aa]];

        output->molecular_type_tag[aa] = molecular_type_tag[order_of_final_arrays[aa]];

        output->atom_ID[aa] = atom_ID_tag[order_of_final_arrays[aa]];

        output->bonded_to[2*aa] = bonded_to[2*order_of_final_arrays[aa]];
        output->bonded_to[2*aa+1] = bonded_to[2*order_of_final_arrays[aa]+1];

        output->angle[3*aa] = angle[3*aa];
        output->angle[3*aa+1] = angle[3*aa+1];
        output->angle[3*aa+2] = angle[3*aa+2];

        //printf("\natom_ID[%i] = %i", aa, atom_ID[aa]);

        output->par_pos[3*aa] = positions[3*order_of_final_arrays[aa]];
        output->par_pos[3*aa+1] = positions[3*order_of_final_arrays[aa]+1];
        output->par_pos[3*aa+2] = positions[3*order_of_final_arrays[aa]+2];
    }

    for(aa=0; aa<9; aa++)
    {
        output->latt_vec[aa] = input->latt_vec[aa];
    }

    output->scale_factor = input->scale_factor;

    output->num_types = par_list->master_num_types_atoms;

    for(aa=0; aa<par_list->master_num_types_atoms; aa++)
    {
        output->par_num[aa] = par_list->master_actual_num_each_type[aa];
    }

    output->par_count = total_num_atoms;

    for(aa=0; aa<par_list->master_num_types_atoms; aa++)
    {
        copyString(output->par_ident[aa], par_list->master_ident[aa]);
        //printf("\noutput->par_ident[%i] = %s", aa, output->par_ident[aa]);
    }


    //printf("\n\n\n\n");

    free(order_of_final_arrays);
    free(pp_ID);
    free(POSCAR_ID);
    free(atom_ID);
    free(pp_pos);
    free(pp_pos_frac);
    free(positions);
    free(molecular_ID);
}


void write_POSCAR(POSCAR_struct* output, char* file_path, placed_particles* par_list, int finished)
{
    int a, b, c, d;
    FILE* out = NULL;

    int holder;

    if(finished==0)
    {
        out = fopen(file_path, "w");
    }

    int* molecular_type_tag_array = NULL;

    //Beginning tasks needed for writing the resubmittable POSCAR file

    if(finished==1)
    {
        int e, f, g, h, i;
        f=0;

        int resubmit;

        int num_each_type_molecule[par_list->num_types_molecules];

        for(e=0; e<par_list->num_types_molecules; e++)
        {
            num_each_type_molecule[e] = par_list->molecule_info[e].actual_number_placed;
        }


        int max = find_max_int(output->molecular_ID, output->par_count);
        par_list->molecule_count = max+1;


        char* POSCAR_file_path = (char*) malloc(1000*sizeof(char));
        char* final_file_path = (char*) malloc(1000*sizeof(char));
        FILE* master_input_reread = NULL;
        char* master_input_file_path = (char*) malloc(200*sizeof(char));

        sprintf(master_input_file_path, "%smaster_input.txt", file_prefix);

        master_input_reread = fopen(master_input_file_path, "r");

        if(master_input_reread==NULL)
        {
            printf("\nFailed to open master_input for a reread!\n");
            exit(666);
        }
        fscanf(master_input_reread, "Falling particles: %i", &holder);
        fscanf(master_input_reread, "\nResubmit: %i", &resubmit); printf("\nresubmit = %i", resubmit);
        fscanf(master_input_reread, "\nInput POSCAR name: %s", POSCAR_file_path);

        sprintf(final_file_path, "%sRESUBMIT_%s", file_prefix, POSCAR_file_path);

        out = fopen(final_file_path, "w");

        int sum_of_previous = 0;




        molecular_type_tag_array = (int*) malloc(par_list->molecule_count*sizeof(int));



        f=0;
        for(g=0; g<par_list->num_types_molecules; g++)
        {
            for(h=0; h<par_list->molecule_info[g].actual_number_placed; h++)
            {
                molecular_type_tag_array[f] = g;
                f++;
            }
        }
    }


    ///////////////////////////////////////////////////////////////////////////////////


    if(out==NULL)
    {
        printf("\nFailed to create output POSCAR!");
        exit(2452);
    }

    for(a=0; a<output->num_types; a++)
    {
        fprintf(out,"%s   ", output->par_ident[a]); //Printing header information
    }

    fprintf(out,"\n");

    fprintf(out, "    %.16lf", output->scale_factor); //writing the scale factor

    for(a=0; a<3; a++) //writing lattice vectors
    {
        fprintf(out,"\n      %.16lf    %.16lf    %.16lf", output->latt_vec[3*a], output->latt_vec[3*a+1], output->latt_vec[3*a+2]);
    }

    fprintf(out,"\n"); //adding a space for formatting purposes

    for(b=0; b<output->num_types; b++) //writing atom identities
    {
        fprintf(out, "   %s", par_list->master_ident[b]);
        //printf("\n\t%s", par_list->master_ident[b]);
    }

    fprintf(out, "\n ");

    for(c=0; c<output->num_types; c++) //writing atom number line
    {
        fprintf(out,"   %i", output->par_num[c]);
    }

    fprintf(out, "\nDirect");

    for(d=0; d<output->par_count; d++)
    {
        if(finished==0)
        {
            fprintf(out, "\n  %.16lf    %.16lf    %.16lf    #%i %i", output->par_pos[3*d], output->par_pos[3*d+1], output->par_pos[3*d+2], output->molecular_ID[d], output->molecular_type_tag[d]);
        }
        else if((finished==1) && (output->molecular_ID[d]==-1))
        {
            fprintf(out, "\n  %.16lf    %.16lf    %.16lf    #%i %i", output->par_pos[3*d], output->par_pos[3*d+1], output->par_pos[3*d+2], output->molecular_ID[d], -1);
        }
        else if((finished==1) && (output->molecular_ID[d]>-1))
        {
            fprintf(out, "\n  %.16lf    %.16lf    %.16lf    #%i %i", output->par_pos[3*d], output->par_pos[3*d+1], output->par_pos[3*d+2], output->molecular_ID[d], molecular_type_tag_array[output->molecular_ID[d]]);
        }

    }


    fclose(out);

}


int return_radius_index(char idents[3])
{
    int radius_index;
    int a;


    for(a=0; a<num_elements; a++)
    {
        if(strcmp(element_ID[a], idents)==0)
        {
            radius_index = a;
            a = num_elements;
        }
    }

    //printf("\nelement_ID[%i] = %s\tcovalent = %lf\tvdw = %lf", radius_index, element_ID[radius_index], covalent_radii[radius_index], vdw_radii[radius_index]);

    return radius_index;
}


void write_XDATCAR(char* file_prefix, int num_POSCARS, placed_particles* master, int iteration)
{
    /*

    1) Use read_poscar() to collect position information from each previously written POSCAR
    2) Create initial lines, position lines (don't forget Konfig=(12 spaces)#)
    3) Close first POSCAR
    4) Open next POSCAR
    5) Skip to position lines and copy them over
    6) Repeat
    7) ???
    8) Profit

    */

    /// num_POSCARS is actually the number of POSCARS after POSCAR 0.
    /// The total number of POSCARS is num_POSCARS+1

    int a, b;

    FILE* xdatcar = NULL;

    int length = 0;

        length = strlen(file_prefix) + strlen("POSCAR_0.XDATCAR") +1;

        char* xdatcar_file_path = (char*) malloc(length*sizeof(char));

        sprintf(xdatcar_file_path, "%s%s", file_prefix, "POSCAR_0.XDATCAR");

    if(iteration==0)
    {
        xdatcar = fopen(xdatcar_file_path, "w");
    }
    else if(iteration>0)
    {
        xdatcar = fopen(xdatcar_file_path, "a");
    }


    if(xdatcar==NULL)
    {
        printf("\nFailed to create XDATCAR file!");
        exit(666);
    }

    int current_POSCAR = iteration;

    POSCAR_struct* holder = (POSCAR_struct*) malloc(sizeof(POSCAR_struct));


    char* POSCAR_file_path = (char*) malloc( (strlen(file_prefix)+strlen("POSCAR_")+strlen(".POSCAR")+sizeof(current_POSCAR))*sizeof(char) );





    sprintf(POSCAR_file_path, "%s%s%i%s", file_prefix, "POSCAR_", current_POSCAR, ".POSCAR");

    holder->num_types = master->master_num_types_atoms;

    read_poscar(holder, POSCAR_file_path, 0); //gathering information for header lines of XDATCAR

    if(iteration==0)
    {
        for(a=0; a<master->master_num_types_atoms; a++)
        {
            fprintf(xdatcar,"%s   ", holder->par_ident[a]); //printing element ID line
        }

        fprintf(xdatcar,"\n\t%0.16lf", holder->scale_factor);

        for(a=0; a<3; a++)
        {
            fprintf(xdatcar, "\n\t  %0.16lf\t%0.16lf\t%0.16lf", holder->latt_vec[3*a], holder->latt_vec[3*a+1], holder->latt_vec[3*a+2]);
        }
    }

    //start loop here

    fprintf(xdatcar,"\nKonfig=            %i", current_POSCAR+1);

    for(b=0; b<holder->par_count; b++)
    {
        fprintf(xdatcar, "\n  %0.16lf\t%0.16lf\t%0.16lf\t#", holder->par_pos[3*b], holder->par_pos[3*b+1], holder->par_pos[3*b+2]);
    }


    free(POSCAR_file_path);
    free(holder);
    fclose(xdatcar);
}


void write_LAMMPS(POSCAR_struct* input, placed_particles* par_list, int num_original_POSCAR_atom_types, char file_extension[])
{
     //TASK 1: Calculate number of different atom types
     int num_atoms = input->par_count; //filled
     int num_bonds = 0; //filled
     int num_angles = 0; //filled
     int num_dihedrals = 0; //filled
     int num_impropers = 0; /**Ask about impropers-- setting to 0 until I find out what they are**/
     int num_atom_types; //filled
     int num_bond_types = par_list->num_types_molecules; //for first pass, set equal to num_types_molecules (assuming one type per molecule)
     int num_angle_types = par_list->num_types_molecules; //for first pass, set equal to num_types_molecules (assuming one type per molecule)

     /**I'll need to update num_bond_types and num_angle_types with actual values that are determined from each molecules' structure**/

     int num_dihedral_types = 0; //for first pass, set equal to 0 (assuming water is the only molecule used)
     int num_improper_types = 0; /**Ask about impropers-- setting to 0 until I find out what they are**/


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


     double cos_alpha;
     double cos_beta;
     double cos_gamma;

     double lx, ly, lz, xy, xz, yz;

     double a_length = input->scale_factor*dist(0,0,0, input->latt_vec[0], input->latt_vec[1], input->latt_vec[2]); printf("\na_length = %lf\n\n\n\n", a_length);
     double b_length = input->scale_factor*dist(0,0,0, input->latt_vec[3], input->latt_vec[4], input->latt_vec[5]);
     double c_length = input->scale_factor*dist(0,0,0, input->latt_vec[6], input->latt_vec[7], input->latt_vec[8]);

     double x_vector[3] = {input->latt_vec[0], input->latt_vec[1], input->latt_vec[2]};
     double y_vector[3] = {input->latt_vec[3], input->latt_vec[4], input->latt_vec[5]};
     double z_vector[3] = {input->latt_vec[6], input->latt_vec[7], input->latt_vec[8]};


     unitize(x_vector);
     unitize(y_vector);
     unitize(z_vector);


    int* bond_amounts = (int*) malloc(num_bond_types*sizeof(int));
    int* angle_amounts = (int*) malloc(num_angle_types*sizeof(int));


    //Calculating num_atoms

     int a, b, c, d, e, f, g;

     for(a=0; a<par_list->num_types_molecules; a++)
     {
         bond_amounts[a] = (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms-1);
         angle_amounts[a] = (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms-2);

         num_bonds += (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms-1);
         num_angles += (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms-2);
         num_dihedrals += (par_list->molecule_info[a].actual_number_placed)*(par_list->molecule_info[a].structure->num_atoms-3);
     }

     if(num_bonds<0)
     {
         num_bonds = 0;
     }

     if(num_angles<0)
     {
         num_angles = 0;
     }

     if(num_dihedrals<0)
     {
         num_dihedrals = 0;
     }


     //Calculating num_atoom_types
     ///Check back on this section later

    num_atom_types = num_original_POSCAR_atom_types;

    for(b=0; b<par_list->num_types_molecules; b++)
    {
        num_atom_types+=par_list->molecule_info[b].num_atom_types;
    }




    //Calculating geometric stuff for unit cell

    cos_alpha = dot(y_vector, z_vector);
    cos_beta = dot(x_vector, z_vector);
    cos_gamma = dot(x_vector, y_vector);

    lx = a_length;
    xy = b_length*cos_gamma;
    xz = c_length*cos_beta;
    ly = sqrt((b_length*b_length) - (xy*xy));
    yz = ((b_length*c_length*cos_alpha)-(xy*xz))/ly;
    lz = sqrt(c_length*c_length-xz*xz-yz*yz);

    //Creating "Atoms" section

    int* POSCAR_type_ident = (int*) malloc(input->par_count);
    int current = 0;

    int first_molecule_type = input->num_types;
    int last_molecule_type = par_list->master_num_types_atoms - input->num_types;

    for(d=0; d<input->num_types; d++)
    {
        for(e=0; e<input->par_num[d]; e++)
        {
            POSCAR_type_ident[current] = d+1;
        }
    }

    double* cartesian_coords = frac_to_cartesian(input->latt_vec, input->par_pos, input->par_count, input->scale_factor);

    char** atom_types = (char**) malloc(num_atoms*sizeof(char*));

    char** master_atom_types = (char**) malloc(num_atom_types*sizeof(char*));

    int* atom_type_number = (int*) malloc(num_atoms*sizeof(int));



    //Printing "Atoms" section

int ident_index = 0;
int current_num = input->par_num[ident_index];
int current_type = 0;

    for(c=0; c<num_atoms; c++)
    {
        if(c==current_num)
        {
            ident_index++;
            current_num+=input->par_num[ident_index];
            //printf("\nINCREMENTING @ c = %i\n", c);
        }

        atom_types[c] = (char*) malloc(20*sizeof(char));

        sprintf(atom_types[c], "%s_%i",input->par_ident[ident_index], input->molecular_type_tag[c]);

        if(c==0)
        {
            master_atom_types[current_type] = (char*) malloc(20*sizeof(char));
            copyString(master_atom_types[current_type], atom_types[c]);
        }
        else
        {
            if(strcmp(master_atom_types[current_type], atom_types[c])!=0)
            {
                current_type++;
                master_atom_types[current_type] = (char*) malloc(20*sizeof(char));
                copyString(master_atom_types[current_type], atom_types[c]);
            }
        }
    }

    for(c=0; c<num_atoms; c++)
    {
        for(d=0; d<num_atom_types; d++)
        {
            if(strcmp(master_atom_types[d], atom_types[c])==0)
            {
                atom_type_number[c] = d+1;
            }
        }
    }

    FILE* filePtr = NULL;

    char file_name[128];

    sprintf(file_name, "%sdata%s", file_prefix, file_extension);

    filePtr = fopen(file_name, "w");

    fprintf(filePtr, "LAMMPS data file. CGCMM style. atom_style full generated by MCPliQ (the quantum wrapper)");
    fprintf(filePtr, "\n %i atoms", num_atoms);
    fprintf(filePtr, "\n %i bonds", num_bonds);
    fprintf(filePtr, "\n %i angles", num_angles);
    fprintf(filePtr, "\n %i dihedrals", num_dihedrals);
    fprintf(filePtr, "\n %i impropers", num_impropers);
    fprintf(filePtr, "\n %i atom types", num_atom_types);
    fprintf(filePtr, "\n %i bond types", num_bond_types);
    fprintf(filePtr, "\n %i angle types", num_angle_types);
    fprintf(filePtr, "\n %i dihedral types", num_dihedral_types);
    fprintf(filePtr, "\n %i improper types", num_improper_types);

    fprintf(filePtr, "\n 0.000000 %lf xlo xhi", lx);
    fprintf(filePtr, "\n 0.000000 %lf ylo yhi", ly);
    fprintf(filePtr, "\n 0.000000 %lf zlo zhi", lz);
    fprintf(filePtr, "\n %lf %lf %lf xy xz yz\n\n", xy, xz, yz);

    fprintf(filePtr, " Pair Coeffs\n");

    for(c=0; c<num_atom_types; c++)
    {
        fprintf(filePtr, "\n %i X Y # %s", c+1, master_atom_types[c]);
    }

    fprintf(filePtr, "\n\n Bond Coeffs\n");

    for(c=0; c<num_bond_types; c++)
    {
        fprintf(filePtr, "\n %i X Y # molecule %i", c+1, c);
    }

    fprintf(filePtr, "\n\n Angle Coeffs\n");

    for(c=0; c<num_angle_types; c++)
    {
        fprintf(filePtr, "\n %i X Y # molecule %i", c+1, c);
    }

    fprintf(filePtr, "\n\n Masses\n");

    for(c=0; c<num_atom_types; c++)
    {
        fprintf(filePtr, "\n %i X # %s", c+1, master_atom_types[c]);
    }

    fprintf(filePtr, "\n\n Atoms\n");


    for(c=0; c<num_atoms; c++)
    {
        fprintf(filePtr, "\n%i 1 %i X %lf %lf %lf # %s", c+1, atom_type_number[c], input->par_pos[3*c], input->par_pos[3*c+1], input->par_pos[3*c+2], atom_types[c]);
    }


    fprintf(filePtr, "\n\n Bonds\n");

    int current_pos = 0;

    int bond_counter = 0;
    int angle_counter = 0;

    int bonds_written = bond_amounts[bond_counter];
    int angles_written = angle_amounts[angle_counter];


    for(c=0; c<num_atoms; c++)
    {
        if(current_pos>=bonds_written)
        {
            bond_counter++;
            bonds_written+=bond_amounts[bond_counter];
        }

        if(input->bonded_to[2*c]!=input->bonded_to[2*c+1])
        {
            fprintf(filePtr, "\n%i %i %i %i", current_pos+1, bond_counter+1, input->bonded_to[2*c], input->bonded_to[2*c+1]);
            current_pos++;
        }
    }


    fprintf(filePtr, "\n\n Angles\n");

    current_pos = 0;

    for(c=0; c<num_atoms; c++)
    {
        if(current_pos>=angles_written)
        {
            angle_counter++;
            angles_written+=angle_amounts[angle_counter];
        }

        if(input->angle[3*c]!=-1)
        {
            fprintf(filePtr, "\n%i %i %i %i %i", current_pos+1, angle_counter+1, input->angle[3*c], input->angle[3*c+1], input->angle[3*c+2]);
            current_pos++;
        }
    }

    fprintf(filePtr, "\n\n");

    fclose(filePtr);
}




  //////////////////////////////////////////////
 ////    MC looping structure functions    ////
//////////////////////////////////////////////

void energy_calc_simple(POSCAR_struct* surface, int iteration_num, int falling_particles)
{
double energy = 0;

    if(falling_particles==0)
    {
        char* command = (char*) malloc(500*sizeof(char));

        sprintf(command, "cp POSCAR_%i.POSCAR POSCAR", iteration_num);
        printf("\n\n\n\n\n\t\t\t\tPOSCAR_%i.POSCAR (the number is the iteration number)\n\n\n\n\n", iteration_num);


        system(command);


        printf("\nbefore VASP call\tITERATION %i-----------------------------------------------------------\n\n", iteration_num);
        system("mpiexec -n 16 /common/curium/SOFTWARE/executables/vaspph7a"); //calls VASP

        system("rm -f core*");
        printf("\nafter VASP call\t ITERATION %i-----------------------------------------------------------------", iteration_num);

        //This is for Linux systems (Palmetto)

        FILE *system_output = NULL;

        system_output = popen("grep sigma OUTCAR | tail -1 | awk '{print $7}'", "r"); //greps energy from latest run

        if(system_output==NULL)
        {
            printf("\n\n\n\n\n\n\n\n\n\nFailed to read system energy with GREP command. Exiting.\n\n\n\n\n\n\n\n\n\n");
            exit(111);
        }

        //Read the output a line at a time - output it

        fscanf(system_output, "%lf", &energy);
        printf("\n\n\n\n\n\t\t\tENERGY = %lf\n\n\n\n\n\n", energy);

    //    while (fgets(energy, sizeof(energy)-1, system_output) != NULL)
    //    {
    //        printf("%lf", energy);
    //    }

        //close
        pclose(system_output);
    }
    else if(falling_particles==1)
    {
        int a;

        for(a=0; a<surface->par_count; a++)
        {
            energy+=surface->par_pos[3*a+2];
        }
    }

    surface->energy = energy;

    //printf("\nsystem energy = %lf gigatons of TNT\n", surface->energy);
}


POSCAR_struct* random_move(POSCAR_struct* first_output, POSCAR_struct* second_output, placed_particles* replacement, placed_particles* particle_list, double prob_of_BAP, double prob_of_rota, double prob_of_trans, double min_scale, double max_scale, double z_min, double z_max, int iteration_num, int zone, double maximum_translation, int rm_termination_limit)
{
    int total_molecules = 0;
    int rm_termination_counter = 0;
    double* number = num_gen(1, particle_list);
    int* old_pos_array_indices = NULL;
    POSCAR_struct* written = (POSCAR_struct*) malloc(sizeof(POSCAR_struct));
        copy_POSCAR(written, second_output); //copies second_output to written
        written->rm_terminated = 0; //sets initial value of the random move termination counter to 0
    int selection; //selection will be the same as the molecular_ID tag of the molecule. It's used to select a fluid molecule to move randomly
    int intramolecular_selection; //used to determine the molecule type of the selected molecule
    int num_atoms_per_molec; //the number of atoms per molecule in the molecule type of the selected molecule
    int molecule_type; //the molecule type, as determined with intramolecular_selection, of the selected molecule
    int last_molecule_molecular_ID = 0; //stores the molecular_ID of the last molecule of the selected molecule's molecule type
    int a, b, c, d;
    int current_index = 0;
    int last_current_index = 0;

    for(a=0; a<particle_list->num_types_molecules; a++)
    {
        total_molecules+=particle_list->molecule_info[a].actual_number_placed;
    }

    selection = *number*total_molecules; //randomly selects a fluid molecule to replace. Uses index numbers (0->inf)

    intramolecular_selection = selection+1; //used for determining which molecule type has been selected. Uses counting numbers (1->inf)

    for(b=0; b<particle_list->num_types_molecules; b++) //determines which molecular type has been selected
    {
        last_molecule_molecular_ID+=particle_list->molecule_info[b].actual_number_placed;

        if(intramolecular_selection>particle_list->molecule_info[b].actual_number_placed)
        {
            intramolecular_selection-=particle_list->molecule_info[b].actual_number_placed;
        }
        else //this stops the loop after the selected fluid molecule's molecule type has been found, before all molecule types have been run through (in case the molecule type is not the last type in the list)
        {
            molecule_type = b;
            b = particle_list->num_types_molecules; //breaks out of loop once condition is satisfied
        }
    }

    last_molecule_molecular_ID--; //determines
    //printf("\nlast_molecule_molecular_ID = %i", last_molecule_molecular_ID);

    replacement->rand_seed = particle_list->rand_seed; //passes the random seed between structures

    num_atoms_per_molec = particle_list->molecule_info[molecule_type].structure->num_atoms;

    //printf("\nnum_atoms_per_molec = %i\n", num_atoms_per_molec);

    old_pos_array_indices = (int*) malloc(num_atoms_per_molec*sizeof(int));

    int* last_pos_array_indices = (int*) malloc(num_atoms_per_molec*sizeof(int));

    int* new_pos_array_indices = (int*) malloc(num_atoms_per_molec*sizeof(int));


    for(c=0; c<first_output->par_count; c++) //this loop determines the array positions of the atoms in the selected molecule
    { //printf("\nfirst_output->molecule_ID[%i] = %i", c, first_output->molecular_ID[c]);
        if(first_output->molecular_ID[c]==selection)
        {
            old_pos_array_indices[current_index] = c;
            //printf("\nold_pos[%i] = %lf", 3*old_pos_array_indices[current_index], first_output->par_pos[3*old_pos_array_indices[current_index]]);
            //printf("\nold_pos[%i] = %lf", 3*old_pos_array_indices[current_index]+1, first_output->par_pos[3*old_pos_array_indices[current_index]+1]);
            //printf("\nold_pos[%i] = %lf", 3*old_pos_array_indices[current_index]+2, first_output->par_pos[3*old_pos_array_indices[current_index]+2]);
            //printf("\nold_pos_array_indices[%i] = %i", current_index, old_pos_array_indices[current_index]);
            current_index++;
        }
/*
        if(first_output->molecular_ID[c]==last_molecule_molecular_ID) //this if-statement does literally nothing. last_pos_array_indices and new_pos_array_indices aren't used at all later in this function or in the rest of the program
        {
            last_pos_array_indices[last_current_index] = c; //are these two lines even used for anything after this loop? I don't think they are. They are not. A+++++++++++ for extra effort.
            new_pos_array_indices[last_current_index] = c+particle_list->molecule_info[molecule_type].structure->num_atoms;
            //printf("\nlast_pos_array_indices[%i] = %i", last_current_index, last_pos_array_indices[last_current_index]);
            //printf("\nnew_pos_array_indices[%i] = %i", last_current_index, new_pos_array_indices[last_current_index]);
            last_current_index++;
        }
*/
    }

////////////////////////////////////////////////

    number = num_gen(1,particle_list);

//Filling the passed placed_particles structure
    replacement->num_types_molecules = 1;


//Creating and filling the associated molecule_coords structure

    replacement->molecule_info = (molecule_coords*) malloc(sizeof(molecule_coords));

    replacement->molecule_info->max_number = 1;
    replacement->molecule_info->structure = particle_list->molecule_info[molecule_type].structure;
    replacement->molecule_info->num_atom_types = particle_list->molecule_info[molecule_type].num_atom_types;
    replacement->molecule_info->num_each_type = (int*) malloc((replacement->molecule_info->num_atom_types)*sizeof(int));
    replacement->molecule_info->max_num_each_type = (int*) malloc((replacement->molecule_info->num_atom_types)*sizeof(int));
    replacement->molecule_info->replacing = 1;
    replacement->molecule_info->replaced_molecular_ID = selection+last_molecule_molecular_ID;
    replacement->molecule_info->rotate_about_this_atom = *number*num_atoms_per_molec; //I created a duplicate variable (possibly) in particle_placement_shell()

    for(d=0; d<replacement->molecule_info->num_atom_types; d++)
    {
        replacement->molecule_info->num_each_type[d] = particle_list->molecule_info[molecule_type].num_each_type[d];
        replacement->molecule_info->max_num_each_type[d] = particle_list->molecule_info[molecule_type].num_each_type[d];
        copyString(replacement->molecule_info->type_ident[d], particle_list->molecule_info[molecule_type].type_ident[d]);
    }

    replacement->molecule_info->replaced_par_pos_indices = old_pos_array_indices;

///////////////////////////////////////////////////////////////



//Calling functions needed to replace the removed molecule


    if((*number*100)<prob_of_BAP)
    {
        //Bond angle perturbation
        replacement->molecule_info->translate = 0;
        replacement->molecule_info->rotate = 0;
        replacement->molecule_info->BAP = 1;
    }
    else if((*number*100)<prob_of_BAP+prob_of_rota)
    {
        //Rotation
        replacement->molecule_info->translate = 0;
        replacement->molecule_info->rotate = 1;
        replacement->molecule_info->BAP = 0;
    }
    else if((*number*100)<=prob_of_BAP+prob_of_rota+prob_of_trans)
    {
        //Translation
        replacement->molecule_info->translate = 1;
        replacement->molecule_info->rotate = 0;
        replacement->molecule_info->BAP = 0;
    }
    else
    {
        printf("\nERROR: invalid move probabilties in random_move()\nExiting");
        exit(4743);
    }



//////////////////

    //printf("\nfile_path = %s", file_path);

    organize_atom_idents(replacement, second_output); //This fills in the master arrays in the particle_list* structure

    written->rm_terminated = particle_placement_shell(replacement, second_output, zone, min_scale, max_scale, z_min, z_max, maximum_translation); //This is the first attempt at a random move

    rm_termination_counter+=written->rm_terminated; //This counts the number of times a random move has been attempted unsuccessfully

    while((written->rm_terminated==1) && (rm_termination_counter<rm_termination_limit)) //this loop continues attempting a random move until it's successful
    {
        written->rm_terminated = particle_placement_shell(replacement, second_output, zone, min_scale, max_scale, z_min, z_max, maximum_translation); printf("\n4162"); //This continues attempting a random move
        rm_termination_counter+=written->rm_terminated; //printf("\nrm_termination_counter = %i", rm_termination_counter); //This counts the number of times a random move has been attempted unsuccessfully
    }

    //printf("\nin random_move(): x1 = %lf  y1 = %lf  z1 = %lf", replacement->molecule_info->coords[0], replacement->molecule_info->coords[1], replacement->molecule_info->coords[2]);
    //printf("\nin random_move(): x2 = %lf  y2 = %lf  z2 = %lf", replacement->molecule_info->coords[3], replacement->molecule_info->coords[4], replacement->molecule_info->coords[5]);
    //printf("\nin random_move(): x3 = %lf  y3 = %lf  z3 = %lf", replacement->molecule_info->coords[6], replacement->molecule_info->coords[7], replacement->molecule_info->coords[8]);


    if(replacement->terminate_replacing==0) //this checks to see if the random move was successful at all
    {
        replacement_sorter(replacement, second_output, written, last_molecule_molecular_ID, selection); //If the random move was successful, then it sends off to replacement_sorter() to incorporate the new molecule into the main POSCAR_struct
    }

//////////////////






    free(number);
    free(old_pos_array_indices);
    //free(last_pos_array_indices);
    //free(new_pos_array_indices);




    return written;
}


void replacement_sorter(placed_particles* replacement, POSCAR_struct* second_output, POSCAR_struct* written, int last_molecular_ID, int selection)
{

    /*

    Steps:

    1) Find the molecular_ID indices (and calculate the position array indices) of the atoms in the molecule selected for replacement and the last_molecular_ID atoms

    2) Create a replaced_atom_ID array. This will be used to record the element of each of the replaced atoms in the order they're found in the position/molecular_ID arrays,
          by referencing the replacement structure's master_ident and recording the master_ident index for each type encountered.
          To determine which type is encountered, use the molecule_coords structure's type_ident, num_atom_types, and num_each_type.

    3) Create a new_atom_ID array by referencing replacement's master_ident. The order of this array will be the order of atom placement/the internal coordinates file atom order.
          This will be used to record the element of each of the new atoms.

    4) Using replaced_atom_ID, new_atom_ID, and a new array called matches, complete the following procedure:
        a) Cycle through replaced_atom_ID in the outer for() loop
        b) Cycle through new_atom_ID in the second for() loop
        c) Compare the current replaced_atom_ID to the current new_atom_ID
        d) If a match is found, record the current new_atom_ID index in the current matches[] position
        e) Increment the current matches[] index by 1
        f) Once each new_atom_ID has been matched to a replaced_atom_ID of the same element, move on to the next step

    5) Now, we have every index for new_atom_ID (as many atoms as there are in the specific molecule being replaced) matched to its corresponding atom in the replaced molecule
          Using this, we can

    */

    int* replaced_molecular_ID_pos = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int)); //filled
    int* last_molecular_ID_pos = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int)); //filled

    int* replaced_atom_ID = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int)); //filled
    int* new_atom_ID = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int)); //filled

    int* matches = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int));

            //basically, matches is the order of placing the new atom in the old molecular_ID and position arrays
            //matches[0] contains the replaced_molecular_ID_pos of the atom the first atom in the new molecule will replace.
            //Ex: with H20, oxygen is placed first. so, if the replaced oxygen's molecular_ID array index is 77, then
            //   its replaced_molecular_ID_pos value (replaced_molecular_ID_pos[0]) will be 77,
            //   so matches[0] will be 77, meaning that the first oxygen in the new molecule will be replacing the 78th (77+1, adjusting for index->counting numbers) atom in the parent POSCAR


    int a, b, c, d, e, f;
    int replaced_molecule_index = 0;
    int last_molecular_ID_pos_index = 0;
    int new_atom_ID_index = 0;

    int num_atoms_already_counted = 0;

    int* chosen_already = (int*) malloc(replacement->molecule_info->structure->num_atoms*sizeof(int));

    for(a=0; a<replacement->molecule_info->structure->num_atoms; a++)
    {
        chosen_already[a] = 0;
    }

    for(a=0; a<second_output->par_count; a++) //this loop completes step 1:
    {
        if(second_output->molecular_ID[a]==selection)
        {
            replaced_molecular_ID_pos[replaced_molecule_index] = a; //records molecular_ID[] array position of each replaced atom

            for(b=0; b<second_output->num_types; b++)
            {
                if((a>=num_atoms_already_counted) && (a<num_atoms_already_counted+second_output->par_num[b]))
                {
                    replaced_atom_ID[replaced_molecule_index] = b; //this completes step 2
                }
                num_atoms_already_counted+=second_output->par_num[b];
            }
            num_atoms_already_counted = 0;
            replaced_molecule_index++;
        }

        if(second_output->molecular_ID[a]==last_molecular_ID)
        {
            last_molecular_ID_pos[last_molecular_ID_pos_index] = a; //records the molecular_ID[] array position of each atom of the last molecule
            //printf("\nlast_molecular_ID_pos[%i] = %i", last_molecular_ID_pos_index, last_molecular_ID_pos[last_molecular_ID_pos_index]);
            last_molecular_ID_pos_index++;
        }
    }


    for(b=0; b<replacement->molecule_info->structure->num_atoms; b++) //filling new_atom_ID (completing step 3)
    {
        for(c=0; c<replacement->master_num_types_atoms; c++)
        {
            if(strcmp(replacement->master_ident[c], replacement->molecule_info->structure->atom_ident[b])==0)
            {
                new_atom_ID[new_atom_ID_index] = c;
                //printf("\nnew_atom_ID[%i] = %i", new_atom_ID_index, new_atom_ID[new_atom_ID_index]);
                new_atom_ID_index++;
            }
        }
    }

//comparing new_atom_ID and replaced_atom_ID


    for(d=0; d<replacement->molecule_info->structure->num_atoms; d++) //loops through new molecule
    {
        for(e=0; e<replacement->molecule_info->structure->num_atoms; e++) //loops through replaced molecule
        {
            if(new_atom_ID[d]==replaced_atom_ID[e])
            {
                if(chosen_already[e]==1)
                {
                    //printf("\nduplicate!");
                    e++;
                }
                chosen_already[e] = 1;
                matches[d] = replaced_molecular_ID_pos[e];
                //printf("\nd = %i\te = %i\treplaced_molecular_ID_pos[%i] = %i",d, e, e, replaced_molecular_ID_pos[e]);
                e = replacement->molecule_info->structure->num_atoms; //breaks out of loop to prevent overwriting of first match with the last match found
                //Now I have to prevent the loop from stopping on the first instance of two or more of the same molecule type /
                /**^^^^^^^^^^^^^^^^^^^^^^^^^^???????????????????????????????/**/
            }
        }
    }

//printf("\nmatches[2] = %i", matches[2]);
//position array indices for replaced atom = x) 3*matches[0]
                                        //   y) 3*matches[0]+1
                                        //   z) 3*matches[0]+2


    double* fractional_coords = (double*) malloc(3*replacement->molecule_info->structure->num_atoms*sizeof(double));
//printf("\n\nx_c = %lf\ty_c = %lf\tz_c = %lf", replacement->molecule_info->coords[0], replacement->molecule_info->coords[1], replacement->molecule_info->coords[2]);


    fractional_coords = cartesian_to_frac(second_output->latt_vec, replacement->molecule_info->coords, replacement->molecule_info->structure->num_atoms, second_output->scale_factor);
//printf("\nnum_atoms_per_molec = %i", replacement->molecule_info->structure->num_atoms);
//printf("\nin replacement_sorter():");
//printf("\nx_f = %lf\ty_f = %lf\tz_f = %lf", fractional_coords[0], fractional_coords[1], fractional_coords[2]);
   //printf("\nx2 = %lf  y2 = %lf  z2 = %lf", fractional_coords[3], fractional_coords[4], fractional_coords[5]);
   //printf("\nx3 = %lf  y3 = %lf  z3 = %lf", fractional_coords[6], fractional_coords[7], fractional_coords[8]);

    for(f=0; f<replacement->molecule_info->structure->num_atoms; f++)
    {
        if(fractional_coords[3*f]>1)
        {
            fractional_coords[3*f]-=1;
        }
        else if(fractional_coords[3*f]<0)
        {
            fractional_coords[3*f]+=1;
        }

        if(fractional_coords[3*f+1]>1)
        {
            fractional_coords[3*f+1]-=1;
        }
        else if(fractional_coords[3*f+1]<0)
        {
            fractional_coords[3*f+1]+=1;
        }

        if(fractional_coords[3*f+2]>1)
        {
            fractional_coords[3*f+2]-=1;
        }
        else if(fractional_coords[3*f+2]<0)
        {
            fractional_coords[3*f+2]+=1;
        }
        written->par_pos[3*matches[f]] = fractional_coords[3*f];
        written->par_pos[3*matches[f]+1] = fractional_coords[3*f+1];
        written->par_pos[3*matches[f]+2] = fractional_coords[3*f+2];
        //printf("\nmatches[%i] = %i", f, matches[f]);
        //printf("\nmatches[%i] = %i  x%i = %lf  y%i = %lf  z%i = %lf",f, matches[f], f+1, fractional_coords[3*f], f+1, fractional_coords[3*f+1], f+1, fractional_coords[3*f+2]);
        //printf("\nmatches[%i] = %i  x%i = %lf  y%i = %lf  z%i = %lf",f, matches[f], f+1, written->par_pos[3*matches[f]], f+1, written->par_pos[3*matches[f]+1], f+1, written->par_pos[3*matches[f]+2]);
    }


    free(replaced_molecular_ID_pos);
    free(last_molecular_ID_pos);
    free(replaced_atom_ID);
    free(new_atom_ID);
    //free(matches);
    free(chosen_already);
    //free(fractional_coords);

}



  ////////////////////////////////
 ////    Helper functions    ////
////////////////////////////////


void copyString(char *to, char* from)
{
    for( ; *from!='\0'; ++from, ++to)
        *to = *from;

        *to= '\0';
}


double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double distance;

    distance = pow(( ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)) + ((z1-z2)*(z1-z2)) ),0.5);

    return distance;
}


double* num_gen(int runs, placed_particles* par_list)
{

    //Function description--------------------------------------------------------------
    /*
        -This function is used to generate pseudorandom numbers between 0 and 1.

        -This function accepts the number of numbers needed

        -This function returns pseudorandom numbers in an array of length = runs.
    */


    double* x=(double *)calloc(runs+50, sizeof(double));
    double* test=(double *)calloc(2, sizeof(double));
    double* vals=(double *)calloc(runs,sizeof(double));
    int m=1;

    int a;

    int first = 0;

    double rand_seed = par_list->rand_seed;

    //printf("\nruns = %i\n",runs);

    //Sets up output file--------------------------------------------------------------
/*
    FILE *values;
    values=fopen("D:\\values.csv","w");

    if(values==NULL)
    {
        printf("\nFailed to open file for values pointer!\n");
    }
*/


    //Checks allocation------------------------------------------------------------------
    if(x==NULL)
    {
        printf("\n\nFailed to allocate memory for pointer variable x!\n\n");
        exit(1112);
    }

    if(test==NULL)
    {
        printf("\n\nFailed to allocate memory for pointer variable test!\n\n");
        exit(1118);
    }

    if(vals==NULL)
    {
        printf("\nFailed to allocate memory for vals!\n");
        exit(1124);
    }


    //Generates pseudo-random numbers---------------------------------------------------

    if(rand_seed==initial_rand_seed)
    {
        first = 1;
    }


    while(m<(runs+50))
    {

        if((first==0) && (m==1))
        {
            test[1] = rand_seed;
            x[1] = rand_seed;
            m = 2;

        }


        if ((m==1))
        {
            test[0]=(double)cos(100*time(NULL));
            if((test[0]<=1)&&(test[0]>=0))
            {
                x[m]=test[0];
                test[1]=test[0];
                //printf("\n\tx[%i] = %f",m,x[m]);
                //fprintf(values,"%lf\n",x[m]);
            }
            else
            {
                m--;
            }
        }
        else if ((m>1)&&(pow(x[m-1]*x[m-1],0.5)<=0.5))
        {
            test[0]=((((((double)(pow(pow(-test[1]*-test[1],0.5),0.5))*((200000)*sin((time(NULL)*-test[1]*pow(eul,test[1]*test[0]))))+(10000))/40000))+0.26));
            if((test[0]<=1)&&(test[0]>=0))
            {
                x[m]=test[0];
                test[1]=test[0];
                //printf("\n\tx[%i] = %f",m,x[m]);
                //fprintf(values,"%lf\n",x[m]);
            }
            else
            {
                m--;
            }
        }
        else if ((m>1)&&(pow(x[m-1]*x[m-1],0.5)>0.5))
        {
            test[0]=-((((((double)(pow(pow(test[1]*test[1],0.5),0.5))*((200000)*sin((time(NULL)*-test[1]*pow(eul,test[1]*-test[0]))))+(10000))/40000))+0.26));
            if((test[0]<=1)&&(test[0]>=0))
            {
                x[m]=test[0];
                test[1]=test[0];
                //printf("\n\tx[%i] = %f",m,x[m]);
                //fprintf(values,"%lf\n",x[m]);
            }
            else
            {
                m--;
            }
        }

        m++;

    }

    //Places all values from x (after the 50th number) into vals for return to calling function----------------------
    for(a=0; a<runs; a++)
    {
        vals[a]=x[a+50];
        //printf("\nvals[%i] = %f",a,vals[a]);
    }

    par_list->rand_seed = vals[runs-1];

    //Closes file, frees memory, ends function---------------------------------------------------------------------

    //fclose(values);

    free(test);
    free(x);

    return vals;
}


double dot(double* vector1, double* vector2) //accepts two 3-component arrays (vectors)
{
    double dotProd = (vector1[0]*vector2[0]) + (vector1[1]*vector2[1]) + (vector1[2]*vector2[2]);

    return dotProd;
}


double* create_unit_vector(double* head_pos, double* tail_pos) //accepts two 3-component arrays (positions)
{
    double* unit_vec = (double*) malloc(3*sizeof(double));
    double length = dist(head_pos[0], head_pos[1], head_pos[2], tail_pos[0], tail_pos[1], tail_pos[2]);

    unit_vec[0] = (head_pos[0] - tail_pos[0]) / length;
    unit_vec[1] = (head_pos[1] - tail_pos[1]) / length;
    unit_vec[2] = (head_pos[2] - tail_pos[2]) / length;

    return unit_vec;
}


int isEqual(int x, double y) //this is used to compare a floating point value to an integer value while accounting for rounding errors due to the storage methods used
{
    double allowable_error = 0.00001;

//printf("\n\tx = %i\ty = %lf", x, y);
    if(((x-y)<allowable_error) && ((x-y)>(-1*allowable_error)))
    {//printf("\n\tEqual");
        return 1;
    }
    else
    {//printf("\n\tNot equal");
        return 0;
    }
}


void copy_POSCAR(POSCAR_struct* dest, POSCAR_struct* source)
{
    //Send an allocated instance of POSCAR_struct for dest (i.e., make the declaration POSCAR_struct* dest = (POSCAR_struct*) malloc(sizeof(POSCAR_struct)) before calling this function)

    dest->par_num = (int*) malloc(source->num_types*sizeof(int));
    dest->par_pos = (double*) malloc(3*source->par_count*sizeof(double));
    dest->molecular_ID = (int*) malloc(source->par_count*sizeof(int));
    dest->molecular_type_tag = (int*) malloc(source->par_count*sizeof(int));
    dest->atom_ID = (int*) malloc(source->par_count*sizeof(int));
    dest->bonded_to = (int*) malloc(2*source->par_count*sizeof(int));
    dest->angle = (int*) malloc(3*source->par_count*sizeof(int));

    int a, b, c;

    dest->scale_factor = source->scale_factor;

    for(a=0; a<9; a++)
    {
        dest->latt_vec[a] = source->latt_vec[a];
    }

    dest->num_types = source->num_types;

    for(b=0; b<dest->num_types; b++)
    {
        dest->par_num[b] = source->par_num[b];
        copyString(dest->par_ident[b], source->par_ident[b]);
    }

    dest->par_count = source->par_count;

    for(c=0; c<dest->par_count; c++)
    {
        dest->par_pos[3*c] = source->par_pos[3*c];
        dest->par_pos[3*c+1] = source->par_pos[3*c+1];
        dest->par_pos[3*c+2] = source->par_pos[3*c+2];

        dest->molecular_ID[c] = source->molecular_ID[c];
        dest->molecular_type_tag[c] = source->molecular_type_tag[c];
        dest->atom_ID[c] = source->atom_ID[c];

        dest->bonded_to[2*c] = source->bonded_to[2*c];
        dest->bonded_to[2*c+1] = source->bonded_to[2*c+1];

        dest->angle[3*c] = source->angle[3*c];
        dest->angle[3*c+1] = source->angle[3*c+1];
        dest->angle[3*c+2] = source->angle[3*c+2];
    }

    dest->energy = source->energy;

    dest->rm_terminated = source->rm_terminated;

}


double* cross(double a_pos[3], double b_pos[3])
{
    double* product = (double*) malloc(3*sizeof(double));

    product[0] = b_pos[2]*a_pos[1] - b_pos[1]*a_pos[2];         // i
    product[1] = -1*(b_pos[2]*a_pos[0] - b_pos[0]*a_pos[2]);    // j
    product[2] = b_pos[1]*a_pos[0] - b_pos[0]*a_pos[1];

    return product;
}


void unitize(double vector[3])
{
    double length = dist(0,0,0, vector[0], vector[1], vector[2]);

    vector[0] = vector[0]/length;
    vector[1] = vector[1]/length;
    vector[2] = vector[2]/length;
}


void VASP_energies(POSCAR_struct* blah, int iteration_num, int move_type, int accepted, int accepted_move_number, int attempted_move_number, int total_num_translate, int total_num_rotate, int finished, double probability_of_rotation)
{
    FILE* energies = NULL;
    FILE* excel_energies = NULL;

    char file_path[1000];
    char excel_file_path[1000];

    sprintf(file_path,"%senergies.txt", file_prefix);
    sprintf(excel_file_path, "%sexcel_energies.txt", file_prefix);

    if(iteration_num==0)
    {
        energies = fopen(file_path, "w");
        excel_energies = fopen(excel_file_path, "w");
        fprintf(energies, "Init = Initial configuration\t\tT = Translation\t\tR = Rotation\t\tBAP = Bond Angle Perturbation\n\n");
    }
    else
    {
        energies = fopen(file_path, "a");
        excel_energies = fopen(excel_file_path, "a");
    }

    if(energies==NULL)
    {
        printf("\n\n3) Failed to create file energies.txt\n\n");
    }

    if(excel_energies==NULL)
    {
        printf("\n\n4) Failed to create file excel_energies.txt\n\n");
    }

    //fprintf(energies, "\nEnergy = %lf\tIteration = %i", blah->energy, iteration_num);



    char* accepted_string = malloc(4*sizeof(char));

    if(accepted==1)
    {
        accepted_string = "yes";
    }
    else if(accepted==0)
    {
        accepted_string = "no";
    }


    char* move_string = malloc(15*sizeof(char));

    if(move_type==0)
    {
        move_string = "T";
    }
    else if(move_type==1)
    {
        move_string = "R";
    }

    if(iteration_num==0)
    {
        move_string = "Init";
    }

    double translations = (double) total_num_translate;
    double rotations = (double) total_num_rotate;

    double alpha, beta, percent_T, percent_R;


    if(finished==0)
    {
        if(accepted==1)
        {
            fprintf(excel_energies, "%i\t\t%lf\n", accepted_move_number, blah->energy);
        }

        if((iteration_num+1)%10==0)
        {
            fprintf(energies, "%i\t\tEnergy: %lf\t\tAccepted?: %s\t\tAccepted Moves: %i\t\tMove type: %s\t\tTranslations: %i\tRotations: %i\n\n", iteration_num+1, blah->energy, accepted_string, accepted_move_number, move_string, total_num_translate, total_num_rotate);
        }
        else
        {
            fprintf(energies, "%i\t\tEnergy: %lf\t\tAccepted?: %s\t\tAccepted Moves: %i\t\tMove type: %s\t\tTranslations: %i\tRotations: %i\n", iteration_num+1, blah->energy, accepted_string, accepted_move_number, move_string, total_num_translate, total_num_rotate);
        }
    }
    else if(finished==1)
    {
        percent_T = (100*translations)/(rotations+translations);
        percent_R = (100*rotations)/(rotations+translations);

        fprintf(energies, "\n\n\n\nACCEPTED MOVE STATISTICS\n\nPercent R = %lf\tPercent T = %lf\tPercent BAP = %lf\n\n\n\n", percent_R, percent_T, 0);
    }

    fclose(energies);
    fclose(excel_energies);
}


int find_max_int(int array[], int num_members)
{ //I can make this a general purpose maximum finder by having it return the array index of the largest element instead of the largest element itself.
    int max = array[0];
    int a;

    for(a=1; a<num_members; a++)
    {
        if(array[a]>max)
        {
            max = array[a];
        }
    }

    return max;
}


  ///////////////////////////////////////////////
 ////    Coordinate Conversion Functions    ////
///////////////////////////////////////////////

double* frac_to_cartesian(double latt_vec[], double pos_frac[], int p_count, double scale_factor)
{
    //Function description---------------------------------------------------------
    /*
        -This function is used to convert frac. lattice coordinates to cartesian coordinates

        -This function accepts frac. lattice coordinates

        -This function returns cartesian coordinates
    */
/*
    FILE* fileptr = NULL;

    fileptr = fopen("D://frac_to_cartesian.txt", "w");

    if(fileptr==NULL)
    {
        printf("\nFailed to create frac_to_cartesian.txt!");
        exit(88);
    }
*/

    int a;

    double* p_carte = NULL;
    p_carte = (double*) malloc(3*p_count*sizeof(double));

    if(p_carte==NULL)
    {
        printf("\np_count = %i", p_count);
        printf("\nFailed to create p_carte!\n");
        exit(8);
    }


    for(a=1; a<=p_count; a++)
    {
        p_carte[3*a-3]=scale_factor*((pos_frac[3*a-3]*latt_vec[0])+(pos_frac[3*a-2]*latt_vec[3])+(pos_frac[3*a-1]*latt_vec[6]));
        p_carte[3*a-2]=scale_factor*((pos_frac[3*a-3]*latt_vec[1])+(pos_frac[3*a-2]*latt_vec[4])+(pos_frac[3*a-1]*latt_vec[7]));
        p_carte[3*a-1]=scale_factor*((pos_frac[3*a-3]*latt_vec[2])+(pos_frac[3*a-2]*latt_vec[5])+(pos_frac[3*a-1]*latt_vec[8]));

        //fprintf(fileptr,"%0.16lf   %0.16lf   %0.16lf\n", p_carte[3*a-3], p_carte[3*a-2], p_carte[3*a-1]);
  //      printf("\nx_f = %lf\ty_f = %lf\t z_f = %lf", pos_frac[3*a-3], pos_frac[3*a-2], pos_frac[3*a-1]);
    }

    //fclose(fileptr);
    return p_carte;
}

double* cartesian_to_frac(double latt_vec[], double pos_cart[], int p_count, double scale_factor)
{//Unit cell transformation

    double* p_frac = NULL;
    p_frac = (double*) malloc(3*p_count*sizeof(double));

    if(p_frac==NULL)
    {
        printf("\nFailed to create p_frac!\n");
        exit(9);
    }

    double adj[3][3]; // convention used: rows x cols
    double det = 0;
    double inverse[3][3];

    double a11, a12, a13, a21, a22, a23, a31, a32, a33;

    int a, b, c;



    //Calculating the matrix adjugate and determinant
/*    printf("\nscale_factor = %lf", scale_factor);
    for(b=0; b<9; b++)
    {
        printf("\nlatt_vec[b] = %lf", latt_vec[b]);
    }
*/

    a11 = scale_factor*latt_vec[0]; //defining the M matrix
    a12 = scale_factor*latt_vec[1];
    a13 = scale_factor*latt_vec[2];

    a21 = scale_factor*latt_vec[3];
    a22 = scale_factor*latt_vec[4];
    a23 = scale_factor*latt_vec[5];

    a31 = scale_factor*latt_vec[6];
    a32 = scale_factor*latt_vec[7];
    a33 = scale_factor*latt_vec[8];

    adj[0][0] = (a22*a33)-(a23*a32); //defining the adjugate
    adj[0][1] = (a13*a32)-(a12*a33);
    adj[0][2] = (a12*a23)-(a13*a22);

    adj[1][0] = (a23*a31)-(a21*a33);
    adj[1][1] = (a11*a33)-(a13*a31);
    adj[1][2] = (a11*a23)-(a31*a21);

    adj[2][0] = (a21*a32)-(a22*a31);
    adj[2][1] = (a12*a31)-(a11*a32);
    adj[2][2] = (a11*a22)-(a12*a21);


    det = (a11*a22*a33)+(a12*a23*a31)+(a12*a21*a32)-(a12*a21*a33)-(a11*a23*a32)-(a13*a22*a31);
//printf("\n\n\n\ndeterminant = %lf\n", det);

    inverse[0][0] = (1/det)*adj[0][0]; //calculating the inverse of M (also called P)
    inverse[0][1] = (1/det)*adj[1][0];
    inverse[0][2] = (1/det)*adj[2][0];

    inverse[1][0] = (1/det)*adj[0][1];
    inverse[1][1] = (1/det)*adj[1][1];
    inverse[1][2] = (1/det)*adj[2][1];

    inverse[2][0] = (1/det)*adj[0][2];
    inverse[2][1] = (1/det)*adj[1][2];
    inverse[2][2] = (1/det)*adj[2][2];

/*
printf("\n\ninverse[0][0] = %lf", inverse[0][0]);
printf("\ninverse[0][1] = %lf", inverse[0][1]);
printf("\ninverse[0][2] = %lf", inverse[0][2]);

printf("\n\ninverse[1][0] = %lf", inverse[1][0]);
printf("\ninverse[1[1] = %lf", inverse[1][1]);
printf("\ninverse[1][2] = %lf", inverse[1][2]);

printf("\n\ninverse[2][0] = %lf", inverse[2][0]);
printf("\ninverse[2][1] = %lf", inverse[2][1]);
printf("\ninverse[2][2] = %lf", inverse[2][2]);
*/

    for(a=1; a<=p_count; a++)
    {
        p_frac[3*a-3] = ((inverse[0][0]*pos_cart[3*a-3]) + (inverse[0][1]*pos_cart[3*a-2]) + (inverse[0][2]*pos_cart[3*a-1]));
        p_frac[3*a-2] = ((inverse[1][0]*pos_cart[3*a-3]) + (inverse[1][1]*pos_cart[3*a-2]) + (inverse[1][2]*pos_cart[3*a-1]));
        p_frac[3*a-1] = ((inverse[2][0]*pos_cart[3*a-3]) + (inverse[2][1]*pos_cart[3*a-2]) + (inverse[2][2]*pos_cart[3*a-1]));
    }

    return p_frac;
}



  /////////////////////////////////////////
 ////    Memory Deleting Functions    ////
/////////////////////////////////////////


void free_placed_particles(placed_particles* par_list)
{
    int a;

    for(a=0; a<par_list->num_types_molecules; a++)
    {   //freeing members of the molecule_struct
        free(par_list->molecule_info[a].structure->bonded_to);
        free(par_list->molecule_info[a].structure->lengths_bond);
        free(par_list->molecule_info[a].structure->third_atom);
        free(par_list->molecule_info[a].structure->angles_bond);
        free(par_list->molecule_info[a].structure->fourth_atom);
        free(par_list->molecule_info[a].structure->dihedral);
        free(par_list->molecule_info[a].structure);
        //freeing members of the molecule_info structure(s)
        free(par_list->molecule_info[a].min_distance_array);
        free(par_list->molecule_info[a].num_each_type);
        free(par_list->molecule_info[a].max_num_each_type);
        free(par_list->molecule_info[a].coords);
        free(par_list->molecule_info[a].molecular_ID);
        free(par_list->molecule_info[a].replaced_par_pos_indices);
    }

    free(par_list->molecule_info);
    free(par_list->master_actual_num_each_type);
    free(par_list);
}



void free_POSCAR_struct(POSCAR_struct* poscar_file)
{
    free(poscar_file->par_num);
    free(poscar_file->par_pos);
    free(poscar_file->molecular_ID);
    free(poscar_file);
}







