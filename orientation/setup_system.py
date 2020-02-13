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

    ## Below is for testing ##
    # Lambda = U.T.dot(I.dot(U))
    #
    # print(Lambda)
    # print(np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0))

    #print("")
    #print(U)

    return U

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
