from orientation.setup_system import get_universe, read_stride

import MDAnalysis as mda
import os
import pickle
import filecmp

gro_file = os.path.join(os.getcwd(), "data", "b3_syst_protein_only.gro")
traj_file = os.path.join(os.getcwd(), "data", "b3_frm_human_b1_r0_400ns_noPBCWhole_noJump_Center_SKIP10.xtc")

stride_file = os.path.join(os.getcwd(), "data", "beta3_stride_file.txt")
sec_struc_file = os.path.join(os.getcwd(), "data", "sec_struc.txt")
protein_length = 166

def test_get_universe():

    uni = get_universe(gro_file, traj_file)

    assert isinstance(uni, mda.Universe)
    

def test_get_universe_notraj():

    uni = get_universe(gro_file)
    
    assert isinstance(uni, mda.Universe)

def test_read_stride():

    protein_dictionary = read_stride(stride_file, protein_length, sec_struc_file)

    current_dict = protein_dictionary

    # compare to ref
    ref_dict_path = os.path.join(os.getcwd(), "data", "ref_stride_output_dict.pkl")

    with open(ref_dict_path, 'rb') as file:
        ref_dict = pickle.load(file)

    assert current_dict == ref_dict

