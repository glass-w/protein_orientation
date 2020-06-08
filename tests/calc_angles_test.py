import os

import numpy as np

import MDAnalysis as mda
from orientation.calc_angles import (
    dir_cosine,
    get_com,
    get_principal_axes,
    make_direction_cosine_matrix,
)
from orientation.setup_system import get_universe

sel = "name CA and resid 1:123"

gro_file = os.path.join(os.getcwd(), "data", "b3_syst_protein_only.gro")

traj_file = os.path.join(
    os.getcwd(), "data", "b3_frm_human_b1_r0_400ns_noPBCWhole_noJump_Center_SKIP10.xtc"
)

uni = get_universe(gro_file, traj_file)


def test_get_com():

    ref_com = uni.select_atoms(sel).center_of_mass()

    com = get_com(uni, sel)

    assert np.allclose(ref_com, com)


def test_pa():

    # Methodology taken from Beckstein
    # https://stackoverflow.com/questions/49239475/
    # how-to-use-mdanalysis-to-principal-axes-and-moment-of-inertia-with-a-group-of-at/49268247#49268247

    sel = "name CA and resid 1:123"
    uni = get_universe(gro_file, traj_file)

    CA = uni.select_atoms(sel)
    I = CA.moment_of_inertia()

    U = get_principal_axes(uni, sel)

    Lambda = U.T.dot(I.dot(U))

    assert np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0)


def test_dir_cosine():

    v1 = np.array([0, 0, 1])
    v2 = np.array([0, 1, 0])

    # gives dot prod
    dc = dir_cosine(v1, v2)

    assert dc == 0.0


def test_make_direction_cosine_matrix():

    ref_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    axes = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    CM = make_direction_cosine_matrix(ref_array, axes)

    # convert matrix to array for assertion
    # CM = np.asarray(CM)

    # assert shape and that CM is an identity matrix (for these orthogonal vectors)
    assert (CM.shape[0] == CM.shape[1]) and (CM == np.eye(CM.shape[0])).all()
