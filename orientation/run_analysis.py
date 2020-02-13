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


def write_out(current_pitch_angle, current_roll_angle, angle_list, system, out_file_name, current_frame):
    ''' This is primarily for if the user wants to write out .gro files (i.e. snapshots) if the angles calcualted are between 
    a certain min and max value '''

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
