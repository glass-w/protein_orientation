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

#cython modules
#from dir_angle_new import dir_angle_new as dir_angle

scale_factor = 20 # for drawing the vectors if needed
cut_off = 5 # dor protein - lipid interactions 

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

# not needed anymore?
def dir_angle(v1, v2):

    # Calculate the direction cosine angle (see here for more info: https://en.wikipedia.org/wiki/Direction_cosine)

    direction_cosine_angle = np.arccos(np.dot(np.abs(v1), v2) / np.linalg.norm(v1))

    return direction_cosine_angle


def dir_cosine(v1, v2):

    ''' 
    Given two vectors (v1 and v2) work out the direction cosine between
    them. For more information see: https://en.wikipedia.org/wiki/Direction_cosine)
    '''

    direction_cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

    return direction_cosine


def make_direction_cosine_matrix(ref, axes_set):

    '''
    Given a set of two axes (ref and axes_set) work out the direction cosine matrix
    between them.
    '''

    # unit vectors, take this as the reference set (i.e. the first set of PAs at the first frame)
    ex = ref[0, :]
    ey = ref[1, :]
    ez = ref[2, :]

    # principal axes calcualted at this instance
    u = axes_set[0, :] # 1st PA
    v = axes_set[1, :] # 2nd PA
    w = axes_set[2, :] # 3rd PA

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

# def make_direction_cosine_matrix(axes_set):
#
#     ex = [1, 0, 0]
#     ey = [0, 1, 0]
#     ez = [0, 0, 1]
#
#     u = axes_set[:, 0] # 1st PA
#     v = axes_set[:, 1] # 2nd PA
#     w = axes_set[:, 2] # 3rd PA
#
#     # Work out the angle b/w the 1st PA and the unit vectors
#     u_alpha = dir_angle(u, ex)
#     u_beta = dir_angle(u, ey)
#     u_gamma = dir_angle(u, ez)
#
#     # Work out the angle b/w the 2nd PA and the unit vectors
#     v_alpha = dir_angle(v, ex)
#     v_beta = dir_angle(v, ey)
#     v_gamma = dir_angle(v, ez)
#
#     # Work out the angle b/w the 3rd PA and the unitvectors
#     w_alpha = dir_angle(w, ex)
#     w_beta = dir_angle(w, ey)
#     w_gamma = dir_angle(w, ez)
#
#     # make this into a 3 X 3 matrix and return it
#     c_matrix = np.matrix([[u_alpha, v_alpha, w_alpha],
#                           [u_beta, v_beta, w_beta],
#                           [u_gamma, v_gamma, w_gamma]])
#
#     return c_matrix

def vis_axes(vis, axes_data, center, name):

    # Select the axes based on the input array
    # axis1 = axes_data[:, 0]
    # axis2 = axes_data[:, 1]
    # axis3 = axes_data[:, 2]

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

# cast_check() is not needed anymore using Euler angle formalism?
def cast_check(angle, measurement, vector, pitch_plane, roll_plane):

    '''Takes an angle and works out the sign of that angle in each quadrant'''

    # Check what pitch measurement we want
    if measurement == 'pitch':

        # Check what plane we want this measurement to be in
        if pitch_plane == 'zx' or 'xz':

            opp = vector[2] # the opposite is the z coord of the vector
            adj = vector[0] # the adjacent is the x coord of the vector

        elif pitch_plane == 'zy' or 'yz':

            opp = vector[2]  # the opposite is the z coord of the vector
            adj = vector[1]  # the adjacent is the y coord of the vector


    # Check what roll measurement we want
    elif measurement == 'roll':

        # Check what plane we want this measurement to be in
        if roll_plane == 'zy' or 'yz':

            opp = vector[2] # the opposite is the z coord of the vector
            adj = vector[1] # the opposite is the y coord of the vector

        elif roll_plane == 'zx' or 'xz':

            opp = vector[2]  # the opposite is the z coord of the vector
            adj = vector[1]  # the opposite is the x coord of the vector


    # Calculate the hypotenuse of the triangle defined by the current vector

    hyp = np.sqrt(opp**2 + adj**2)

    # Now calculate the sine, cosine and tangent

    sin_theta = opp / hyp
    cos_theta = adj / hyp
    tan_theta = opp / adj

    # We need the values to be positive, to verify we check them against the CAST diagram (below) and adjust the
    # values accordingly.

    ########## CAST Diagram ##########
    #                                #
    #             90 deg             #
    #               ^                #
    #               |                #
    #          S(2) |  A(1)          #
    #               |                #
    # 180 deg ------ ------> 360 deg #
    #               |                #
    #          T(3) |  C(4)          #
    #               |                #
    #            270 deg             #
    #                                #
    ##################################

    if sin_theta >= 0 and cos_theta <= 0 and tan_theta <= 0: # in second quadrant, 90 -> 180

        # print(measurement + " in 2nd quad: ", "sin(theta) = ", sin_theta, " cos(theta) = ", cos_theta, " tan(theta) = ",
        #       tan_theta, " ANGLE = ", np.rad2deg(angle))

        #angle = np.deg2rad(90) + angle

        angle = np.deg2rad(180) - angle

    elif tan_theta >= 0 and cos_theta <= 0 and sin_theta <= 0: # in third quadrant, 180 -> 270

        # print(measurement + " in 3rd quad: ", "sin(theta) = ", sin_theta, " cos(theta) = ", cos_theta, " tan(theta) = ",
        #       tan_theta, " ANGLE = ", np.rad2deg(angle))

        angle = np.deg2rad(180) + angle

    elif cos_theta >= 0 and tan_theta <= 0 and sin_theta <= 0: # in 4th quadrant, 270 -> 0

        # print(measurement + " in 4th quad: ", "sin(theta) = ", sin_theta, " cos(theta) = ", cos_theta, " tan(theta) = ",
        #       tan_theta, " ANGLE = ", np.rad2deg(angle))

        #angle = np.deg2rad(270) + angle
        angle = np.deg2rad(360) - angle

    return angle


def read_stride(file, num_prots, chain_length, specific_region):
    ''' This reads the output from running stride structure.pdb, this is used to identify beta sheets and alpha
    helices. Due to flexible loops the calculated principal axes can differ so using more stable regions can give
    less noisy data.

    Prior to this function the user must run "stride file.pdb > stride_file.txt"

    '''

    # TODO: work on the specific region option

    x = []
    ss_list = ['E']#, 'G', 'H', 'I']

    ss_list = ['E', 'G'] # strand (E) and 310 helix (G)
    with open(file) as f:

        for line in f.readlines():

            if line.splitlines()[0][0] == 'A':

                if line.splitlines()[0][24] in ss_list:
                    res = (line.splitlines()[0][12], line.splitlines()[0][13], line.splitlines()[0][14])
                    x.append(int(''.join(res)))

    #chain_labels = ['A', 'B', 'C']#, 'D', 'E']

    chain_dict = {'chain ' + str(i): {'resids': [],
                                      'angle_pa1': [],
                                      'angle_pa2': [],
                                      'angle_pa3': [],
                                      } for i in range(num_prots)}

    print (x)
    print(len(x))
    print("HERE")
    print(num_prots, chain_length)

    print(chain_dict)

    for i in range(num_prots):

        if i == 0:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in x if 1 <= t <= chain_length]

        # Need to test below works on a multi chain system
        else:

            chain_dict['chain ' + str(i)]['resids'] = [t for t in x if
                                               (i * chain_length) + 1 <= t <= ((i+1) * chain_length)]  # 124 to 246


    #chain_dict['chain 1']['resids'] = [t for t in x if 1 <= t <= chain_length] # 1 to 123

    #print(chain_dict)

    #chain_dict['chain 2']['resids'] = [t for t in x if (chain_length + 1) <= t <= (2 * chain_length)] # 124 to 246

    #chain_dict['chain 3']['resids'] = [t for t in x if (1 + (2 * chain_length)) <= t <= (3 * chain_length)] # 247 to 369

    print(chain_dict)

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


def run_traj(universe_files, chain_info, skip_val, calc_method, vector_sel_resids, states, run, cut_off, p_plane,
             r_plane, angle_method):

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
        print("Frame = ", ts.frame, ", time = ", ts.time / 1000)
        print("")

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

                # In the first frame we use this basis as the reference for all others
                if i == 0: 
                    ref = pa_array_row_wise

                # Calculate the direction cosine matrix between the current orthonormal basis
                # and the reference. The reference is taken as the basis calcualted in the first frame.

                dir_cosine_matrix = make_direction_cosine_matrix(ref=ref, axes_set=pa_array_row_wise)

                # Create a rotation object from the direction cosine matrix
                rotation_obj = R.from_dcm(dir_cosine_matrix.T)

                #euler_angles1 = rotation_obj1.as_euler('xyz', degrees=True) # extrinsic

                # Intrinisic rotation formalism used
                # Rotate about Z first, then Y, then X. These correspond to Roll, Pitch, and Yaw
                # (https://en.wikipedia.org/wiki/Euler_angles#Alternative_names)
                # i.e. going from user defined basis --> current set of orthonormal vectors at current frame.
                euler_angles_array = rotation_obj.as_euler('ZYX', degrees=True) # intrinsic

                # Store the euler_angles_array so we can keep track of them
                euler_angle_store.append(euler_angles_array)

                #testing
                # print("reference: ", ref)
                # print("directional cosine matrix:", dir_cosine_matrix)
                print("apply rotation to ref: ", rotation_obj.apply(ref))
                print(" ")
                #print("apply rotation to ref: ", rotation_obj1.apply(pa_array_row_wise))
                # print("current basis at frame (row wise): ", pa_array_row_wise)

                print("ref basis ", ref)
                print(" ")
                print("rotated ref ", dir_cosine_matrix @ ref) # '@' used for matrix multiplication
                print(" ")
                print(" basis at frame: ", pa_array_row_wise)
                print(" ")


                print('EULER ANGLES: roll, pitch, yaw')
                print(euler_angles_array)
                print(" ")

                ax1 = pa_array_row_wise[0] # i.e. the roll axis
                ax2 = pa_array_row_wise[1] # i.e. the pitch axis
                ax3 = pa_array_row_wise[2] # i.e. the yaw axis

            # If the user has defined their own axes there is no need to flip them
            # i.e. at every frame three vectors are drawn from the CoM to the three 
            # user defined residues (these may not be / probably won't be orthogonal)
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

    parser.add_argument("-angle_method", dest="angle_method", type=str, default='euler',
                        help="The method you want to use to calculate the angles, euler or plane")

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

    # If the user wants to write out different states (i.e. pitch roll etc) initialise the dicts and angles
    if options.states_list == True:
        calc_states = True

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

    elif options.states_list == False:

        calc_states = False

        states_dict = None









    # what is this doing???
    #######################################################################################################

    # Get the start and end of the protein based on the users selection for the centre of mass
    start_of_protein = int(options.com_selection.split()[-1].split(':')[0])
    end_of_protein = int(options.com_selection.split()[-1].split(':')[1])

    # Load a temp universe, this is to get resids and protein lengths etc to be used in the main loop
    # We assume that the user is analysing multiple repeats of the SAME system
    universe = get_universe(systems[0][0], systems[0][1])

    # select the first universe to make resids
    resid_list = universe.select_atoms("protein").residues.resids

    prot_length = len(
        universe.select_atoms("resid " + str(start_of_protein) + ":" + str(end_of_protein) + " and name CA"))

    chain_dict = read_stride(options.stride_file, options.num_of_proteins, int(prot_length), False)

    #######################################################################################################




    print("Pool Execution")
    start = time.time()

    # results = Parallel(n_jobs=10,
    #                    verbose=3,
    #                    backend="multiprocessing")(
    #     delayed(run_traj)(universe_files=systems,
    #                       chain_info=chain_dict,
    #                       skip_val=options.skip,
    #                       calc_method=options.method,
    #                       vector_sel_resids=options.res_vector_sel,
    #                       states=calc_states,
    #                       run=i,
    #                       cut_off=cut_off,
    #                       p_plane=options.pitch_plane,
    #                       r_plane=options.roll_plane) for i, system in enumerate(systems))


    results = Parallel(n_jobs=options.nprocs,
                       verbose=3,
                       backend="multiprocessing")(
        delayed(run_traj)(universe_files=systems,
                          chain_info=chain_dict,
                          skip_val=options.skip,
                          calc_method=options.method,
                          vector_sel_resids=options.res_vector_sel,
                          states=calc_states,
                          run=i,
                          cut_off=cut_off,
                          p_plane=options.pitch_plane,
                          r_plane=options.roll_plane,
                          angle_method=options.angle_method) for i, system in enumerate(systems))

    delta = time.time() - start
    delta = delta / 60.0
    print("Pool Finished in: " + str(delta))

  
  
  
  
  
   #### BELOW IS ALL TO DO WITH THE STATES CALC ####

    # Get the second element of each item in the results list, this corresponds to the "state" information returned by
    # run_traj
    extracted_state_info = [item[1] for item in results]

    # Populate a new dictionary that now only contains the "state" data
    state_data = {}

    for entry in extracted_state_info:
        state_data.update(entry)

    df_list = []





    if calc_states == True:

       # pitch_ranges = [(0, 90), (90, 180), (180, 270), (270, 360)]

        pitch_ranges = [(-10, 10), (25, 35)]

        min_pitch_list = [x[0] for x in pitch_ranges]

        max_pitch_list = [x[1] for x in pitch_ranges]

        dict_of_ranges = {str(k) : {'frames': [], 'contacts': []} for k in pitch_ranges}

        # Loop over all of the ranges we are interested in
        for i, range_set in enumerate(dict_of_ranges):

            # Loop over every run we have
            for run in state_data:
                # i.e. system 0

                # Loop over the pitch values for the current run
                for j, entry in enumerate(state_data[run]['pitch']):

                    if min_pitch_list[i] <= entry <= max_pitch_list[i]:
                        dict_of_ranges[range_set]['frames'].append(state_data[run]['frame'][j])
                        dict_of_ranges[range_set]['contacts'].append(state_data[run]['contact'][j])

       # print("HELLO")
      #  print(dict_of_ranges)

        # Initiate new dict that will store the averaged interactions for each state
        average_dict_of_ranges = {'state ' + str(pitch_range): [] for pitch_range in pitch_ranges}
        average_dict_of_ranges_2 = {'state ' + str(pitch_range): [] for pitch_range in pitch_ranges}

        for range_set in dict_of_ranges:

            # print("Loop over ranges")
            # print(range_set)

            sum_state = np.zeros((166,))  # need to generalise this
            sum_frame = np.zeros((166,))  # need to generalise this

            for i, contact_array in enumerate(dict_of_ranges[range_set]['contacts']):

                sum_state += np.sum(contact_array, axis=1)

            try:  # since some of the angle ranges may not have been visited and the length of frames will be 0

                sum_state = sum_state / len(dict_of_ranges[range_set]['frames'])

            except ValueError:  # if there is no sampling in this range, skip
                continue

            else:

                average_dict_of_ranges['state ' + str(range_set)].append(sum_state)




            for i, frame in enumerate(dict_of_ranges[range_set]['frames']):

                sum_frame += np.sum(dict_of_ranges[range_set]['contacts'][i], axis=1)

            # try:
            #
            #     sum_frame = sum_state / len(dict_of_ranges[range_set]['frames'])
            #
            # except ValueError:
            #     continue
            #
            # else:

                average_dict_of_ranges_2['state ' + str(range_set)].append([frame, sum_frame])

            # print("Averages")
            # print(average_dict_of_ranges.keys())
            # print(average_dict_of_ranges)

            # sum_state = pd.DataFrame(np.reshape(sum_state, (166, 1)))
            #
            # mask = sum_state.isnull()
            #
            # ax = sns.heatmap(sum_state, mask=mask)
            # ax.set_xlabel('Pitch State ' + str())
            # ax.set_ylabel('Residue Number')
            # plt.show()

           # print("HERE3")
            #print(average_dict_of_ranges)

            df_list.append(pd.DataFrame(average_dict_of_ranges['state ' + str(range_set)]))

           # print(average_dict_of_ranges_2)

        df = pd.concat(df_list)

        # Normalise, take the transpose so that every state is normalised, not every residue
        # Ref: https://stackoverflow.com/questions/41226744/normalize-data-in-pandas-dataframe




        df_norm = (df.T - df.T.min()) / (df.T.max() - df.T.min())

        ax = sns.heatmap(df_norm,  cbar_kws={'label': r'Normalised Average Lipid Distance ($\AA$)'},
                         xticklabels=average_dict_of_ranges.keys(), cmap='magma')

        ax.set_xlabel('Pitch State')
        ax.set_ylabel('Residue Number')
        plt.tight_layout()

        plt.savefig('residue_lipid_interactions_per_pitch_state.svg', dpi=300)
        plt.show()

        # Save the Dataframe in the current directroy

        df.to_pickle('pitch_state_data.pkl')
        df_norm.to_pickle('pitch_state_data_normalised.pkl')

            # for j, pitch_value in enumerate(state_data['system ' + str(i)]):
            #
            #     if min_pitch_list[i] <= pitch_value <= max_pitch_list[i]:
            #
            #         state_dict_range['frames'].append(state_data[run]['frames'])
            #         state_dict_range['frames'].append(state_data[run]['contacts'])