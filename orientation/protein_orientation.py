import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array as dist
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import argparse

from joblib import Parallel, delayed
import time
import sys

import scipy
from scipy.spatial.transform import Rotation as R

# Print Versions
print("MDAnalysis version: ", mda.__version__)
print("NumPy version: ", np.__version__)
print("SciPy version: ", scipy.__version__)
print("Pandas version: ", pd.__version__)

scale_factor = 20 # for drawing the vectors if needed
cut_off = 5 # for protein - lipid interactions 

def get_universe(gro_file, traj_file):

    ''' Load an MDAnalysis universe '''

    u = mda.Universe(gro_file, traj_file)

    return u


def get_principal_axes(universe, selection):

    # Methodology taken from Beckstein
    # https://stackoverflow.com/questions/49239475/
    # how-to-use-mdanalysis-to-principal-axes-and-moment-of-inertia-with-a-group-of-at/49268247#49268247

    # select alpha carbons only
    CA = universe.select_atoms(selection)

    # calculate moment of inertia, and sort eigenvectors
    I = CA.moment_of_inertia()

    # UT = CA.principal_axes()
    # U = UT.T

    values, evecs = np.linalg.eigh(I)
    indices = np.argsort(values)
    U = evecs[:, indices]

    # testing
    # Lambda = U.T.dot(I.dot(U))
    #
    # print(Lambda)
    # print(np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0))

    #print("")
    #print(U)

    return U


def get_com(universe, selection):

    ''' Return the centre of mass of a selection string '''

    com = universe.select_atoms(selection).center_of_mass()

    return com


def dir_cosine(v1, v2):

    ''' 
    Given two vectors (v1 and v2) work out the direction cosine between
    them. For more information see: https://en.wikipedia.org/wiki/Direction_cosine)
    '''

    direction_cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

    return direction_cosine


def make_direction_cosine_matrix(ref, axes_set):

    '''
    Given a set of two axes (ref_option and axes_set) work out the direction cosine matrix
    between them.
    '''
    # Gather the individual vectors of our reference basis
    ex = ref[0, :]
    ey = ref[1, :]
    ez = ref[2, :]

    # principal axes calculated at this instance
    u = axes_set[0, :] # 1st princpal axis (PA)
    v = axes_set[1, :] # 2nd princpal axis (PA)
    w = axes_set[2, :] # 3rd princpal axis (PA)

    # Work out the dot prod b/w the 1st PA and the unit vectors
    u_alpha = dir_cosine(u, ex)
    u_beta = dir_cosine(u, ey)
    u_gamma = dir_cosine(u, ez)

    # Work out the dot prod b/w the 2nd PA and the unit vectors
    v_alpha = dir_cosine(v, ex)
    v_beta = dir_cosine(v, ey)
    v_gamma = dir_cosine(v, ez)

    # Work out the dot prod b/w the 3rd PA and the unit vectors
    w_alpha = dir_cosine(w, ex)
    w_beta = dir_cosine(w, ey)
    w_gamma = dir_cosine(w, ez)

    # make all of the above into a 3 X 3 matrix and return it 
    c_matrix = np.matrix([[u_alpha, u_beta, u_gamma],
                          [v_alpha, v_beta, v_gamma],
                          [w_alpha, w_beta, w_gamma]])                       

    return c_matrix


def vis_axes(vis, axes_data, center, name):

    ''' 
    Visualise the principal axes being used for the calculation.

    vis : the format to write the axes out as (either 'vmd' or 'pymol')
    axes_data :  the array containing the three principal axes
    center : the geometric center of the users selection
    name : the name of the .pdb file (based on the supplied .xtc file)
    '''


    axis1 = axes_data[0]
    axis2 = axes_data[1]
    axis3 = axes_data[2]

    if vis == 'vmd':

        output = open(str(name) + '_pa_vectors.pdb', 'a')

        for i in range(0, (3 * scale_factor)):
            tmp = "ATOM    {0:3d}  CA  ALA A {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(i, i, center[0] + (axis1[0] * i),
                                                                                                    center[1] + (axis1[1] * i),
                                                                                                    center[2] + (axis1[2] * i))
            output.write(tmp)

        output.write("TER\n")

        for j in range(0, (2 * scale_factor)):
            tmp2 = "ATOM    {0:3d}  CA  ALA B {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(j, j, center[0] + (axis2[0] * j),
                                                                                                    center[1] + (axis2[1] * j),
                                                                                                    center[2] + (axis2[2] * j))
            output.write(tmp2)

        output.write("TER\n")


        for k in range(0, (1 * scale_factor)):
            tmp3 = "ATOM    {0:3d}  CA  ALA C {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00\n".format(k, k, center[0] + (axis3[0] * k),
                                                                                                    center[1] + (axis3[1] * k),
                                                                                                    center[2] + (axis3[2] * k))
            output.write(tmp3)

        output.write("TER\nENDMDL\n")

        output.close()

    elif vis == 'pymol':

        # --------------------------------------------------------------------------
        # center axes to the geometric center of the molecule
        # and rescale them by order of eigen values
        # --------------------------------------------------------------------------

        # the large vector is the first principal axis
        point1 = 3 * scale_factor * axis1 + center
        # the medium vector is the second principal axis
        point2 = 2 * scale_factor * axis2 + center
        # the small vector is the third principal axis
        point3 = 1 * scale_factor * axis3 + center

        #pymol_name = pdb_name.replace(".pdb", "_axes.pml")

        pymol_name = ("test_axes.pml")
        with open(pymol_name, "w") as pymol_file:
            pymol_file.write(
                """
                from cgo import *
                axis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                cmd.load_cgo(axis1, 'axis1')
                cmd.load_cgo(axis2, 'axis2')
                cmd.load_cgo(axis3, 'axis3')
                cmd.set('cgo_line_width', 4)
                """ % ( \
                    center[0], center[1], center[2], point1[0], point1[1], point1[2], \
                    center[0], center[1], center[2], point2[0], point2[1], point2[2], \
                    center[0], center[1], center[2], point3[0], point3[1], point3[2]))

        # --------------------------------------------------------------------------
        # create .pml script for nice rendering in Pymol
        # output usage
        # --------------------------------------------------------------------------
        print("\nFirst principal axis (in red)")
        # print("coordinates: ", axis1)
        # print("eigen value: ", eval1)

        print("\nSecond principal axis (in green)")
        # print("coordinates:", axis2)
        # print("eigen value:", eval2)

        print("\nThird principal axis (in blue)")


def read_stride(stride_file, num_prots, chain_length, sec_struc_choice):
    '''
    This reads the output from running stride structure.pdb, this is used to identify beta sheets and alpha
    helices. Due to flexible loops the calculated principal axes can differ so using more stable regions can give
    less noisy data.

    Prior to this function the user must run "stride file.pdb > stride_file.txt" and then use this file for the -stride option.
    '''

    # TODO: work on a specific region option

    sec_struc_dict = {'strand': 'E', '310helix' : 'G', 'alphahelix': 'H'}

    with open(sec_struc_choice, 'r') as ss_file: 
                        choices = [str(choice).rstrip() for choice in ss_file.readline().split(',')] # rstrip() to remove trailing characters

    resid_list = []

    # make a list of the types of secondary structure to include, this is done by referencing sec_struc_dict to get the one letter code as defined by / in the supplied stride file
    sec_struc_list = [sec_struc_dict[sec_feature] for sec_feature in choices]

    print(sec_struc_list)

    # Populate resid_list with the list of residues to use for the main calculation - done with reference to the secondary structure (i.e. ignore flexible loop regions)
    with open(stride_file) as f:
        for line in f.readlines():

            if line.splitlines()[0][0] == 'A': # if at the asignment (ASG) part of the stride file

                if line.splitlines()[0][24] in sec_struc_list: # Read the one letter code column of the stride file
                    res = (line.splitlines()[0][12], line.splitlines()[0][13], line.splitlines()[0][14]) # get the resid (only works for up to 999 currently)
                    resid_list.append(int(''.join(res))) # join them up to make the correct number and append

    # Make the dictionary with the relevant resids and empty lists to store the Euler angles later for each protein in the system
    chain_dict = {'chain ' + str(i): {'resids': [],
                                      'angle_pa1': [],
                                      'angle_pa2': [],
                                      'angle_pa3': [],
                                      } for i in range(num_prots)}


    for i in range(num_prots):

        if i == 0:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in resid_list if 1 <= t <= chain_length]

        # Need to test below works on a multi chain system
        else:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in resid_list if
                                               (i * chain_length) + 1 <= t <= ((i+1) * chain_length)]

    return chain_dict


def write_out(current_pitch_angle, current_roll_angle, angle_list, system, out_file_name, current_frame):
    print("angle list: ", angle_list)

    # Make dict of supplied angle list to check what calculations we want to perform
    angles_to_calc = {'pitch_min': angle_list[0], 'pitch_max': angle_list[1],
                      'roll_min': angle_list[2], 'roll_max': angle_list[3]}

    # Run through the user supplied list of angles and update our dict with the values, converted from a string
    for i, angle in enumerate(angles_to_calc):

        try:
            float(angle_list[i])

        except ValueError:

            angles_to_calc[angle] = None

        else:
            angles_to_calc[angle] = float(angle_list[i])

    prot = system.select_atoms("protein")
    lipid = system.select_atoms("resname POPC")
    prot_lipid = prot + lipid

    #print("ANGLES TO CALC: ", angles_to_calc)

    if (angles_to_calc['roll_min'] and angles_to_calc['roll_max']) is not None: # Only want to analyse roll

        #print("CURRENT ANGLE: ", np.rad2deg(current_roll_angle))

        if angles_to_calc['roll_min'] <= np.rad2deg(current_roll_angle) <= angles_to_calc['roll_max']:

            prot_lipid.write(str(out_file_name) + "_prot_rep_struc_frame" + str(current_frame) + "_roll_range" +
                             ".gro")

    elif (angles_to_calc['pitch_min'] and angles_to_calc['pitch_max']) is not None: # Only want to analyse pitch

        #print("CURRENT ANGLE: ", np.rad2deg(current_pitch_angle))

        if angles_to_calc['pitch_min'] <= np.rad2deg(current_roll_angle) <= angles_to_calc['pitch_max']:

            prot_lipid.write(str(out_file_name) + "_prot_rep_struc_frame" + str(current_frame) + "_pitch_range" +
                             ".gro")

    else: # Write out when both the pitch and roll parameters are satisfied

        print(angles_to_calc)

        if angles_to_calc['pitch_min'] <= np.rad2deg(current_pitch_angle) <= angles_to_calc['pitch_max'] \
                and angles_to_calc['roll_min'] <= np.rad2deg(current_roll_angle) <= angles_to_calc['roll_max']:

            prot_lipid.write(str(out_file_name) + "_prot_rep_struc_frame" + str(current_frame) + "_pitch_and_roll_range"
                                                                                                 + ".gro")


def calc_lipid_contacts(universe, cut_off):

        # reset the centre of mass list for each residue
        com_array = []

        # Cycle through every residue and get it's centre of mass based on it's side chain atoms (from dict)
        for resid in resid_list:
            sel = universe.select_atoms("resid " + str(resid))

            com_array.append(sel.select_atoms("name " + str(term_atom_dict[sel.resnames[0]])).atoms.center_of_mass())

        # convert centre of mass array (a list) for next step, define the head group of the POPC groups.
        com_array = np.array(com_array, dtype="float32")

        popc_head = universe.select_atoms("resname POPC and name P")

        if cut_off:
            # Calc distances of CoMs of residues and POPC heads within a cut off, returns a boolean array
            d = dist(com_array, popc_head.atoms.positions) <= cut_off

        else:
            # Calc distances of CoMs of residues and POPC heads
            d = dist(com_array, popc_head.atoms.positions)

        return d


def run_traj(universe_files, chain_info, skip_val, calc_method, vector_sel_resids, states, run, cut_off, ref_option, ref_basis):

    euler_angle_store = []

    print("Here we go!")

    # find the residues the user specified for the pitch, roll, and yaw
    # these are used to make sure the principal axes vectors calculated 
    # are always pointing in the right direction.
    user_resid_pitch_sel = vector_sel_resids.split(',')[0]
    user_resid_roll_sel = vector_sel_resids.split(',')[1]
    user_resid_yaw_sel = vector_sel_resids.split(',')[2]

    # sort out names for files
    path_name = str(xtc_files[run]).split()[-1].split('.')[-2]

    output_file_name = path_name.split()[-1].split('/')[-1]

    # define the current universe, this is accessed in a for loop. 
    u = mda.Universe(universe_files[run][0], universe_files[run][1])

    # initialise states dictionary - maybe run with if statement?
    states_dict = {'system ' + str(run): {'frame': [], 'pitch': [], 'contact': []}}

    # Go through the current run
    for i, ts in enumerate(u.trajectory[::skip_val]):

        # Show progress
        print("Frame = ", ts.frame, ", Time = ", ts.time / 1000, "ns")

        # At each frame calculate the angle that each chain has moved through
        # here, chain means protein copy... need to change this.

        # TODO: Test on multiple chains, so far have only tested on one chain

        for chain in chain_info:

            # set the resids gathered from the stride file
            # will be used to select the centre of mass, i.e this is the secondary structure of the protein
            sel_resids = ' or resid '.join(map(str, chain_info[chain]['resids']))

            sel = "name CA and (resid " + str(sel_resids) + ")"

            # Define vectors - drawn from centre of mass of body to the resid chosen
            
            # Define the initial pitch vector based on the users choice
            user_ax1 = list(
                u.select_atoms("resid " + str(user_resid_pitch_sel) + " and name CA").atoms.positions[
                    0] - u.select_atoms(sel).center_of_mass())

            # Define the initial roll vector based on the users choice
            user_ax2 = list(
                u.select_atoms("resid " + str(user_resid_roll_sel) + " and name CA").atoms.positions[
                    0] - u.select_atoms(sel).center_of_mass())

            # Define the initial yaw vector based on the users choice
            user_ax3 = list(
                u.select_atoms("resid " + str(user_resid_yaw_sel) + " and name CA").atoms.positions[
                    0] - u.select_atoms(sel).center_of_mass())

            # Normalise the user vectors as they will not be between 0 & 1
            user_ax1 = user_ax1 / np.linalg.norm(user_ax1)
            user_ax2 = user_ax2 / np.linalg.norm(user_ax2)
            user_ax3 = user_ax3 / np.linalg.norm(user_ax3)

            #############################
            #                           # 
            #   Calculate Euler Angles  #  
            #                           #
            #############################

            # if the user has selected a starting set of principal axes to ensure the PAs always point
            # towards them, use these
            if calc_method == "user_pa":

                pa_array = get_principal_axes(u, sel)

                # Make a row wise version for the direction cosine matrix calculation (see the func for
                # why this is needed)
                pa_array_row_wise = pa_array.T

                # Since the way in which the principal axes are calculated in MDAnalysis we need
                # to check that the principal axes at each frame are pointing in the same direction
                # as the one specified by the user, if it is then flip it:

                if np.dot(user_ax1, pa_array_row_wise[0, :]) < 0:

                    pa_array_row_wise[0, :] = pa_array_row_wise[0, :] * -1

                if np.dot(user_ax2, pa_array_row_wise[1, :]) < 0:

                    pa_array_row_wise[1, :] = pa_array_row_wise[1, :] * -1

                if np.dot(user_ax3, pa_array_row_wise[2, :]) < 0:

                    pa_array_row_wise[2, :] = pa_array_row_wise[2, :] * -1

                ##############################
                ###  Get reference basis   ###
                ##############################

                if ref_option == 'first_frame':

                    # In the first frame we use this basis as the reference for all others
                    if i == 0: 
                        ref = pa_array_row_wise

                elif ref_option == 'user':

                    # read each line of the supplied file, corresponds to each basis vector

                    with open(ref_basis, 'r') as file: 
                        a = [int(number) for number in file.readline().split(',')] 
                        b = [int(number) for number in file.readline().split(',')] 
                        c = [int(number) for number in file.readline().split(',')] 

                    ref = np.array([a, b, c])

                else:

                    # use standard box vectors as reference basis
                    ref = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

                ################################
                ###  Calc DCM and Rot matrix ###
                ################################

                # Calculate the direction cosine matrix between the current orthonormal basis and the reference
                dir_cosine_matrix = make_direction_cosine_matrix(ref=ref, axes_set=pa_array_row_wise)

                # Create a rotation object from the direction cosine matrix
                rotation_obj = R.from_dcm(dir_cosine_matrix.T)

                ##############################
                ### Calculate Euler Angles ###
                ##############################

                # Intrinisic rotation formalism used: Rotate about Z first, then Y', then X''. These correspond to Yaw, Pitch, and Roll (https://en.wikipedia.org/wiki/Euler_angles#Alternative_names)
                euler_angles_array = rotation_obj.as_euler('ZYX', degrees=True)

                # Store the euler_angles_array
                euler_angle_store.append(euler_angles_array)

                #testing
                # print("apply rotation to ref: ", rotation_obj.apply(ref))
                # print(" ")
                # print("ref basis ", ref)
                # print(" ")
                # print("rotated ref ", dir_cosine_matrix @ ref) # '@' used for matrix multiplication
                # print(" ")
                # print(" basis at frame: ", pa_array_row_wise)
                # print(" ")

                # Display calculated angles
                yaw, pitch, roll = np.round(euler_angles_array[0], 2), np.round(euler_angles_array[1], 2), np.around(euler_angles_array[2], 2)

                print('EULER ANGLES: Yaw = ' + str(yaw) + ', Pitch = ' + str(pitch) + ', Roll = ' + str(roll) + '\n')

                ax1 = pa_array_row_wise[0] # i.e. the roll axis
                ax2 = pa_array_row_wise[1] # i.e. the pitch axis
                ax3 = pa_array_row_wise[2] # i.e. the yaw axis

            # If the user has defined their own axes there is no need to flip them
            # i.e. at every frame three vectors are drawn from the CoM to the three user defined residues (these may not be / probably won't be orthogonal)
            elif calc_method == 'user':

                ax1 = user_ax1
                ax2 = user_ax2
                ax3 = user_ax3

                    # TODO
                    # calculate angles with Euler again? basis may not be orthogonal...

            # Check if we want to write out the vectors to a pdb file
            if options.vector_traj == True:

                # create coordinates array
                coord = np.array(u.select_atoms(sel).atoms.positions, float)

                # compute geometric center
                center = np.mean(coord, 0)

                vis_axes(vis='vmd', axes_data=[ax1, ax2, ax3], center=center, name=output_file_name)

            #TODO fix this
            if states == True:

                # Calc the contacts
                contacts = calc_lipid_contacts(universe=u, cut_off=cut_off)

                if angle_method == 'plane':

                    # Calc the states

                    states_dict['system ' + str(run)]['frame'].append(ts.frame)
                    states_dict['system ' + str(run)]['pitch'].append(np.rad2deg(pa_1_measured_angle))
                    states_dict['system ' + str(run)]['contact'].append(contacts)

                    # Update pitch and roll dataframes
                    chain_info[chain]['angle_pa1'].append((ts.time / 1000, np.rad2deg(pa_1_measured_angle)))
                    chain_info[chain]['angle_pa2'].append((ts.time / 1000, np.rad2deg(pa_2_measured_angle)))
                    #chain_info[chain]['angle_pa3'].append((ts.time / 1000, np.rad2deg(pa_3_measured_angle)))

                elif angle_method == 'euler':

                    # Calc the states

                    states_dict['system ' + str(run)]['frame'].append(ts.frame)
                    states_dict['system ' + str(run)]['pitch'].append(euler_angles1[1])
                    states_dict['system ' + str(run)]['contact'].append(contacts)

                    # Update pitch and roll dataframes
                    chain_info[chain]['angle_pa1'].append((ts.time / 1000, np.rad2deg(pa_1_measured_angle)))
                    chain_info[chain]['angle_pa2'].append((ts.time / 1000, np.rad2deg(pa_2_measured_angle)))
                    # chain_info[chain]['angle_pa3'].append((ts.time / 1000, np.rad2deg(pa_3_measured_angle)))

        # Save to output files
        #np.savetxt(output_file_name + '_pitch_' + str(p_plane) + '.txt', chain_info[chain]['angle_pa1'])
        #np.savetxt(output_file_name + '_roll_' + str(r_plane) + '.txt', chain_info[chain]['angle_pa2'])

        np.save(output_file_name +'_euler_angles.npy', np.array(euler_angle_store))


    # return chain_info, dict_of_states
    return [chain_info, states_dict]


def init_parser():

    ''' Gather all of the relevant user information '''

    parser = argparse.ArgumentParser(
        description="Calculates the orientation of a user defined region of a protein")

    parser.add_argument("-c", dest="gro_file_list",
                        help='The list of coordinate files [.gro], this takes the form of a text file')

    parser.add_argument("-f", dest="xtc_file_list",
                        help='A list of corrected trajectory files: pbc artifacts removed, no jumping across PBC, this'
                             'is a text file.')

    parser.add_argument("-com_sel", dest="com_selection", type=str,
                        help='The range of resids to use for centre of mass calculation,\
                            in the form of A:B where A and B are integers')

    parser.add_argument("-n", dest="num_of_proteins", type=int,
                        help='The number of protein copies in the system')

    parser.add_argument("-skip", dest="skip", type=int,
                        help="The number of frames to skip", default=1)

    parser.add_argument("-vtraj", dest="vector_traj", type=bool,
                        help="Set to True if you want a trajectory of the vectors", default=False)

    parser.add_argument("-method", dest="method", type=str,
                        help="The vectors can be calculated by 1) a set of user defined vectors based on the centre of\
                            mass of the main selection and the CA of a specified residue OR 2) the method can be used in\
                            combination with (1) and use the principal axes of inertia. In either (1) or (2) the user\
                            must define a set of vectors that roughly correspond to the principal axes - this ensures\
                            that when calculated they always point in the direction specified by the user's vectors\
                            Options: 'user' 'user_pa'")

    parser.add_argument("-res_vector_sel", dest="res_vector_sel", type=str,
                        help="The resids of the residues to use for the roll, pitch, and yaw calculation in the form\
                            A, B, C.")

    parser.add_argument("-stride", dest="stride_file", type=str,
                        help="the name of the stride file to read, a .txt file.\
                            This will be used in combination with the -com_sel selection\
                                to only choose those residues involved in secondary structure.")

    parser.add_argument("-pa_only", dest="pa_single", type=bool,
                        help="If set to True a principal component calculation will be carried out and written to a\
                            .pdb file, this is to help in selecting the appropriate residues for a run.")

    # parser.add_argument("-write_coords", dest="wcoords_list", nargs='+', default=[],
    #                     help="a list of the minimum and maximum pitch and roll angles within which the user wants"
    #                          "to write out the frame. This argument takes a list with four elements in the following"
    #                          "form: min_pitch, max_pitch, min_roll, max_roll. e.g. if you wanted to write out frames"
    #                          "within 0 and 10 degress of pitch and 20 and 30 degress of roll you'd pass the following:"
    #                          "0, 10, 20, 30.")

    # parser.add_argument("-write_states", dest="states_list", type=bool, default=False,
    #                     help="a list of states to write out, here we define pitch and roll as states."
    #                          "We pass a list of minimum and maximum pitch and roll angles within which the user wants")

    parser.add_argument("-nprocs", dest="nprocs", type=int, default=1,
                        help="No. or processes to use, usually set to the number of repeats you have.\
                            Max = max No. of CPUs available to you")

    parser.add_argument("-ref_option", dest="ref_option", type=str, default="standard",
                        help="Choice of what basis of vectors to use as a reference, from which the Euler angles\
                            will be calcualted. Permitted chocies are: 'first_frame', 'user' or 'standard' \
                                where standard is x, y, z = [1,0,0], [0,1,0], [0,0,1].")        

    parser.add_argument("-ref_basis", dest="ref_basis", type=str, default=None,
                        help="The basis vectors to be used as a reference, if not passed the default will be used (see -ref_option).\
                            This should be a .txt file with the x, y, z coordinates one each line e.g.\
                                1,0,0\
                                0,1,0\
                                0,0,1")  

    parser.add_argument("-sec_struc_choice", dest="sec_struc_choice", default=['strand', '310helix', 'alphahelix'],
                        help="A file containing the choice of secondary srtucture to use in the calculation of the centre of mass.\
                            Valid choices include: strand, 310helix, or alphahelix.\
                                In the file these must be comma separated and have no whitespace between them. e.g. strand,310helix")

    return parser.parse_args()

#### MAIN #### 
if __name__ == "__main__":

    # Make dictionary of terminal atoms for each residue,
    # these will be used for the centre of mass calc later
    term_atom_dict = {'ARG': 'NH1 NH2',
                      'HIS': 'ND1 NE2 CE1 CG CD2',
                      'LYS': 'NZ',
                      'ASP': 'OD1 OD2',
                      'GLU': 'OE1 OE2',
                      'SER': 'OG',
                      'THR': 'OG1 CG2',
                      'ASN': 'ND2 OD1',
                      'GLN': 'NE2 OE2',
                      'CYS': 'SG',
                      # 'SEC' : None,
                      'GLY': 'CA',
                      'PRO': 'CA, CB CG CD N',
                      'ALA': 'CB',
                      'VAL': 'CG1 CG2',
                      'ILE': 'CD',
                      'LEU': 'CD1 CD2',
                      'MET': 'CE',
                      'PHE': 'CG, CD2 CE2 CZ CE1 CD1',
                      'TYR': 'OH',
                      'TRP': 'CG CD2 CE3 CZ3 CH2 CZ2 CE2 NE1 CD1'}

    # Get the users options
    options = init_parser()

    # Initialise file lists 
    gro_files = []
    xtc_files = []
    systems = []

    # load .gro files from supplied file
    with open(options.gro_file_list, 'r') as f:
        gro_files = f.read().splitlines()

    # load .xtc files from supplied file
    with open(options.xtc_file_list, 'r') as f:
        xtc_files = f.read().splitlines()


    # populate the systems list with each element a .gro and .xtc file
    # this assumes the list supplied is ordered as sequential repeats. 
    for i in range(len(gro_files)):

        systems.append([gro_files[i], xtc_files[i]])

    # TODO sort this all out, all to do with the states calc...

    # If the user wants to write out different states (i.e. pitch roll etc) initialise the dicts and angles
    # if options.states_list == True:
    #     calc_states = True

        # states_dict = {'system ' + str(i): {'frame': [], 'pitch': [], 'contact': []} for i in range(len(gro_files))}

        # states_dict = {'system ' + str(i): {'frame': [], 'pitch': [], 'contact': []}}

        # Make dict of supplied angle list to check what calculations we want to perform
        # angles_to_calc = {'pitch_min': options.states_list[0], 'pitch_max': options.states_list[1],
        #               'roll_min': options.states_list[2], 'roll_max': options.states_list[3]}
        #
        # # Run through the user supplied list of states and update our dict with the min / max angles, conv from a string
        # for i, angle in enumerate(angles_to_calc):
        #
        #     try:
        #         float(options.states_list[i])
        #
        #     except ValueError:
        #
        #         angles_to_calc[angle] = None
        #
        #     else:
        #         angles_to_calc[angle] = float(options.states_list[i])

    # elif options.states_list == False:

    #     calc_states = False

    #     states_dict = None

    ######################################################
    ### Get information about protein(s) in the system ###
    ######################################################

    # Get the start and end of the protein based on the users selection for the centre of mass
    start_of_protein = int(options.com_selection.split()[-1].split(':')[0])
    end_of_protein = int(options.com_selection.split()[-1].split(':')[1])

    # Load a temp universe, this is to get resids and protein lengths etc to be used in the main loop
    # We assume that the user is analysing multiple repeats of the SAME system
    temp_universe = get_universe(systems[0][0], systems[0][1])
    
    # get resids and length of protein
    resid_list = temp_universe.select_atoms("protein").residues.resids

    prot_length = len(
        temp_universe.select_atoms("resid " + str(start_of_protein) + ":" + str(end_of_protein) + " and name CA"))

    # Initialise dictionary that holds information for each chain in the system (e.g. x3 subunits - N.B. Assumes subunits are all the SAME)
    chain_dict = read_stride(stride_file=options.stride_file, num_prots=options.num_of_proteins, chain_length=int(prot_length), sec_struc_choice=options.sec_struc_choice)

    #####################
    ##                 ##
    ##      RUN        ##    
    ##                 ##
    #####################

    # Run analysis - use multiple processes to run each system at the same time
    print("Pool Execution")
    start = time.time()

    results = Parallel(n_jobs=options.nprocs,
                       verbose=3,
                       backend="multiprocessing")(
        delayed(run_traj)(universe_files=systems,
                          chain_info=chain_dict,
                          skip_val=options.skip,
                          calc_method=options.method,
                          vector_sel_resids=options.res_vector_sel,
                          states=None,
                          run=i,
                          cut_off=cut_off,
                          ref_option=options.ref_option,
                          ref_basis=options.ref_basis) for i, system in enumerate(systems))

    delta = time.time() - start
    delta = delta / 60.0
    print("Pool Finished in: " + str(delta))

  
  
  
  
  

    
