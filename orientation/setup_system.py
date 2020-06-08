import numpy as np

import MDAnalysis as mda


def get_universe(gro_file, traj_file=None):

    """ Load an MDAnalysis universe """

    # print(gro_file, traj_file)

    if traj_file != None:

        u = mda.Universe(gro_file, traj_file)

    else:

        u = mda.Universe(gro_file)

    return u


def read_stride(stride_file, protein_sel_length, sec_struc_choice):
    """
    This reads the output from running stride structure.pdb, this is used to identify beta sheets and alpha
    helices. Due to flexible loops the calculated principal axes can differ so using more stable regions can give
    less noisy data.

    Prior to this function the user must run "stride file.pdb > stride_file.txt" and then use this file for the -stride option.
    """

    # TODO: work on a specific region of the protein as an option

    sec_struc_dict = {"strand": "E", "310helix": "G", "alphahelix": "H"}

    with open(sec_struc_choice, "r") as ss_file:

        choices = [
            str(choice).rstrip() for choice in ss_file.readline().split(",")
        ]  # rstrip() to remove trailing characters

    resid_list = []

    # make a list of the types of secondary structure to include, this is done by referencing sec_struc_dict to get the one letter code as defined by / in the supplied stride file
    sec_struc_list = [sec_struc_dict[sec_feature] for sec_feature in choices]

    # print(sec_struc_list)

    # Populate resid_list with the list of residues to use for the main calculation - done with reference to the secondary structure (i.e. ignore flexible loop regions)
    with open(stride_file) as f:
        for line in f.readlines():

            if (
                line.splitlines()[0][0] == "A"
            ):  # if at the asignment (ASG) part of the stride file

                if (
                    line.splitlines()[0][24] in sec_struc_list
                ):  # Read the one letter code column of the stride file
                    res = (
                        line.splitlines()[0][12],
                        line.splitlines()[0][13],
                        line.splitlines()[0][14],
                    )  # get the resid (only works for up to 999 currently)
                    resid_list.append(
                        int("".join(res))
                    )  # join them up to make the correct number and append

    # Make the dictionary with the relevant resids and empty lists to store the Euler angles later for each protein in the system
    protein_dict = {"resids": [], "angle_pa1": [], "angle_pa2": [], "angle_pa3": []}

    protein_dict["resids"] = [t for t in resid_list if 1 <= t <= protein_sel_length]

    #     # Need to test below works on a multi chain system
    #     else:

    #         chain_dict['chain ' + str(i)]['resids'] = [t for t in resid_list if
    #                                            (i * chain_length) + 1 <= t <= ((i+1) * chain_length)]

    return protein_dict
