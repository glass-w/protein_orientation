from orientation.setup_system import get_universe

import MDAnalysis as mda
import os

gro_file = os.path.join(os.getcwd(), "data", "b3_syst_protein_only.gro")
traj_file = os.path.join(os.getcwd(), "data", "b3_frm_human_b1_r0_400ns_noPBCWhole_noJump_Center_SKIP10.xtc")

def test_get_universe():

    uni = get_universe(gro_file, traj_file)

    assert isinstance(uni, mda.Universe)
    

def test_get_universe_notraj():

    uni = get_universe(gro_file)
    
    assert isinstance(uni, mda.Universe)
