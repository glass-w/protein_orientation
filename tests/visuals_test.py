from orientation.visuals import vis_axes

import numpy as np
import os
import filecmp

axes_array = np.array([[1,0,0],[0,1,0],[0,0,1]])
c = np.array([0,0,0])

def test_vis_axes_vmd():

    vis_axes('vmd', axes_array, c, 'test_vmd')

    ref = os.path.join(os.getcwd(), "data", "test_vmd_pa_vectors.pdb")
    current = os.path.join(os.getcwd(), "", "test_vmd_pa_vectors.pdb")

    assert filecmp.cmp(ref, current)
    
    os.remove(current)
