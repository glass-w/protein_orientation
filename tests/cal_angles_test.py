from setup_system import get_universe
from calc_angles import get_principal_axes

import MDAnalysis as mda
import numpy as np

sel = "name CA and resid 1:123"
uni = get_universe('./data/b3_syst_protein_only.gro')

def test_pa():

    sel = "name CA and resid 1:123"
    uni = get_universe('./data/b3_syst_protein_only.gro')
    
    CA = uni.select_atoms(sel)
    I = CA.moment_of_inertia()

    pa = get_principal_axes(uni, sel)

    Lambda = pa.T.dot(I.dot(U))
    
    assert np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0)
