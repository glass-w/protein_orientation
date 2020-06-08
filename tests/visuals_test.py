import filecmp
import os

import numpy as np

from orientation.visuals import vis_axes

axes_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
c = np.array([0, 0, 0])


def test_vis_axes_vmd():

    vis_axes("vmd", axes_array, c, "test_vmd")

    ref = os.path.join(os.getcwd(), "data", "test_vmd_pa_vectors.pdb")
    current = os.path.join(os.getcwd(), "", "test_vmd_pa_vectors.pdb")

    assert filecmp.cmp(ref, current)

    os.remove(current)


def test_vis_axes_pymol():

    # Currently tested output in Pymol 1.4 only!

    vis_axes("pymol", axes_array, c, "test_pymol")

    ref = os.path.join(os.getcwd(), "data", "test_pymol_pa_vectors.pml")
    current = os.path.join(os.getcwd(), "", "test_pymol_pa_vectors.pml")

    assert filecmp.cmp(ref, current)

    os.remove(current)
