# import modules 
import MDAnalysis as mda
import numpy as np
import scipy
from scipy.spatial.transform import Rotation as R
import sys
import os
import argparse

# Lets run in parallel
from joblib import Parallel, delayed
import time
import sys

# import functions from own modules
from setup_system import get_universe, read_stride
from calc_angles import get_com, make_direction_cosine_matrix, get_principal_axes
from visuals import vis_axes

# Print Versions
print("MDAnalysis version: ", mda.__version__)
print("NumPy version: ", np.__version__)
print("SciPy version: ", scipy.__version__)

def run_single(universe_files, protein_info, calc_method, vector_sel_resids):

    # make visuals of the vectors true, since this is what this function does / is intended for. 
    options.vector_traj == True

    # sort out names for files
    path_name = str(gro_files[0]).split()[-1].split('.')[-2]

    output_file_name = 'REF_' + path_name.split()[-1].split('/')[-1]

    # define the current universe
    u = mda.Universe(universe_files[0])

    # set the resids gathered from the stride file
    # will be used to select the centre of mass, i.e this is the secondary structure of the protein
    sel_resids = ' or resid '.join(map(str, protein_info['resids']))

    sel = "name CA and (resid " + str(sel_resids) + ")"

    # Calculate the raw principal axes - these will be used as a reference for the user
    # as such, we will not need to calculate any angles and just write the vectors out.

    pa_array = get_principal_axes(u, sel)

    # Make a row wise version for the direction cosine matrix calculation (see the func for
    # why this is needed)
    pa_array_row_wise = pa_array.T

    ax1 = pa_array_row_wise[0] # i.e. the roll axis
    ax2 = pa_array_row_wise[1] # i.e. the pitch axis
    ax3 = pa_array_row_wise[2] # i.e. the yaw axis

    ### Write out the PAs ###

    # create coordinates array
    coord = np.array(u.select_atoms(sel).atoms.positions, float)

    # compute geometric center
    center = np.mean(coord, 0)

    vis_axes(vis='vmd', axes_data=[ax1, ax2, ax3], center=center, name=output_file_name)


def run_multiple(universe_files, protein_info, skip_val, calc_method, vector_sel_resids, states, run, ref_option, ref_basis):

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

    # initialise states dictionary - maybe run with if statement? - is this needed anymore?
    states_dict = {'system ' + str(run): {'frame': [], 'pitch': [], 'contact': []}}

    # Go through the current run
    for i, ts in enumerate(u.trajectory[::skip_val]):

        # Show progress
        print("Frame = ", ts.frame, ", Time = ", ts.time / 1000, "ns")

        # set the resids gathered from the stride file
        # will be used to select the centre of mass, i.e this is the secondary structure of the protein
        sel_resids = ' or resid '.join(map(str, protein_info['resids']))

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

        np.save(output_file_name +'_euler_angles.npy', np.array(euler_angle_store))


def init_parser():

    ''' Gather all of the relevant user information '''

    parser = argparse.ArgumentParser(
        description="Calculates the orientation of a user defined region of a protein")

    parser.add_argument("-c", dest="gro_file_list", required=True,
                        help="The list of coordinate files [.gro], this takes the form of a text file with each file location starting on a new line.")

    parser.add_argument("-f", dest="xtc_file_list", required=True,
                        help='The list of corrected trajectory files: pbc artifacts removed, no jumping across PBC. This takes the form of a text file with each file location starting on a new line.')

    parser.add_argument("-com_sel", dest="com_selection", type=str, required=True,
                        help='The range of resids to use for centre of mass calculation, in the form of A:B, where A and B are integers.')

    parser.add_argument("-method", dest="method", type=str, required=True, default="user_pa",
			            help="The vectors can be calculated by 1) a set of user defined vectors based on the centre of mass of the main selection and the alpha carbon (CA) of a specified residue OR 2) the method can be used in combination with (1) and use the principal axes of inertia. In either (1) or (2) the user must define a set of vectors that roughly correspond to the principal axes - this ensures that when calculated they always point in the direction specified by the users vectors. Options: user or user_pa. Default = user_pa")

    parser.add_argument("-n", dest="num_of_proteins", type=int, required=True,
                        help='Number of protein copies in the system, default 1.')

    parser.add_argument("-skip", dest="skip", type=int, default=1,
                        help="The number of frames to skip, default 1.")

    parser.add_argument("-vtraj", dest="vector_traj", type=bool, default=False,
                        help="Set to True if you want a trajectory of the vectors, default False.")

    parser.add_argument("-res_vector_sel", dest="res_vector_sel", type=str, default=None,
                        help="The resids of the residues to use for the roll, pitch, and yaw calculation respectively: in the form A, B, C.")

    parser.add_argument("-stride", dest="stride_file", type=str,
			            help="The name of the stride file to read, a .txt file. This will be used in combination with the -com_sel selection to only choose those residues involved in secondary structure. If using the 'user_pa' method (see -method) this option must be supplied.")

    parser.add_argument("-pa_only", dest="pa_single", type=bool,
                        help="If set to True a principal component calculation will be carried out and written to a .pdb file, this is to help in selecting the appropriate residues for a run. Default False")

    parser.add_argument("-nprocs", dest="nprocs", type=int, default=1,
			            help="Number of processes to use, default=1.")

    parser.add_argument("-ref_option", dest="ref_option", type=str, default="standard",
			            help="Choice of what basis of vectors to use as a reference, from which the Euler angles will be calcualted. Permitted chocies are: 'first_frame', angles will be calculated in reference to the PAs calculated in the first frame. 'user', angles will be calculated in reference to a user defined set of vectors. 'standard' where the standard is x, y, z = [1,0,0], [0,1,0], [0,0,1]. default = 'standard'.")

    parser.add_argument("-sec_struc_choice", dest="sec_struc_choice", default=['strand', '310helix', 'alphahelix'],
			            help="A file containing the choice of secondary structure to use in the calculation of the centre of mass. If using the 'user_pa' method (see -method) this option must be supplied. Valid choices include: 'strand', '310helix', or 'alphahelix'. In the file these must be comma separated and have no whitespace between them. e.g. strand,310helix")

    return parser.parse_args()


if __name__ == "__main__":

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

    ######################################################
    ### Get information about protein(s) in the system ###
    ######################################################

    # Get the start and end of the protein based on the users selection for the centre of mass
    start_of_protein_sel = int(options.com_selection.split()[-1].split(':')[0])
    end_of_protein_sel = int(options.com_selection.split()[-1].split(':')[1])

    # Load a temp universe, this is to get resids and protein lengths etc to be used in the main loop
    # We assume that the user is analysing multiple repeats of the SAME system i.e. multiple MD runs
    temp_universe = get_universe(systems[0][0]) # just the coordinate file
    
    # get resids and length of protein
    resid_list = temp_universe.select_atoms("protein").residues.resids

    prot_sel_length = len(
        temp_universe.select_atoms("resid " + str(start_of_protein) + ":" + str(end_of_protein) + " and name CA"))

    # Initialise dictionary that holds information for protein (resids and angles)
    protein_info = read_stride(stride_file=options.stride_file, protein_sel_length=int(prot_length), sec_struc_choice=options.sec_struc_choice)

    #####################
    ##                 ##
    ##      RUN        ##    
    ##                 ##
    #####################

    if options.pa_single: # If the user just wants to find the principal axes

        single_system = systems[0]

        print("here")
        print(systems[0])
        
        run_single(universe_files=single_system,
                protein_info=protein_info,
                calc_method=options.method,
                vector_sel_resids=options.res_vector_sel)

    else:

        # Run analysis - use multiple processes to run each system at the same time
        print("Pool Execution")
        start = time.time()

        results = Parallel(n_jobs=options.nprocs,
                        verbose=3,
                        backend="multiprocessing")(
            delayed(run_multiple)(universe_files=systems,
                            protein_info=protein_info,
                            skip_val=options.skip,
                            calc_method=options.method,
                            vector_sel_resids=options.res_vector_sel,
                            states=None,
                            run=i,
                            ref_option=options.ref_option,
                            ref_basis=options.ref_basis) for i, system in enumerate(systems))

        delta = time.time() - start
        delta = delta / 60.0
        print("Pool Finished in: " + str(delta))
