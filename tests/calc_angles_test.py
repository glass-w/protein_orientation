from orientation.setup_system import get_universe
from orientation.calc_angles import get_principal_axes

import MDAnalysis as mda
import numpy as np
import os

sel = "name CA and resid 1:123"
gro_file = os.path.join(os.getcwd(), "data", "b3_syst_protein_only.gro")
traj_file = os.path.join(os.getcwd(), "data", "b3_frm_human_b1_r0_400ns_noPBCWhole_noJump_Center_SKIP10.xtc")

uni = get_universe(gro_file, traj_file)

def test_pa():

    sel = "name CA and resid 1:123"
    uni = get_universe('data/b3_syst_protein_only.gro')
    
    CA = uni.select_atoms(sel)
    I = CA.moment_of_inertia()

    pa = get_principal_axes(uni, sel)

    Lambda = pa.T.dot(I.dot(U))
    
    assert np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0)
