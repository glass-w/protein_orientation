import numpy as np


def get_com(universe, selection):

    """ Return the centre of mass of a selection string """

    com = universe.select_atoms(selection).center_of_mass()

    return com


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

    # print("")
    # print(U)

    return U


def dir_cosine(v1, v2):

    """ 
    Given two vectors (v1 and v2) work out the direction cosine between
    them. For more information see: https://en.wikipedia.org/wiki/Direction_cosine)
    """

    direction_cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

    return direction_cosine


def make_direction_cosine_matrix(ref, axes_set):

    """
    Given a set of two axes (ref_option and axes_set) work out the direction cosine matrix
    between them.
    """
    # Gather the individual vectors of our reference basis
    ex = ref[0, :]
    ey = ref[1, :]
    ez = ref[2, :]

    # principal axes calculated at this instance
    u = axes_set[0, :]  # 1st princpal axis (PA)
    v = axes_set[1, :]  # 2nd princpal axis (PA)
    w = axes_set[2, :]  # 3rd princpal axis (PA)

    # Work out the dot prod b/w the 1st PA and the unit vectors
    u_alpha = dir_cosine(u, ex)
    u_beta = dir_cosine(u, ey)
    u_gamma = dir_cosine(u, ez)

    # Work out the dot prod b/w the 2nd PA and the unit vectors
    v_alpha = dir_cosine(v, ex)
    v_beta = dir_cosine(v, ey)
    v_gamma = dir_cosine(v, ez)

    # Work out the dot prod b/w the 3rd PA and the unit vectors
    w_alpha = dir_cosine(w, ex)
    w_beta = dir_cosine(w, ey)
    w_gamma = dir_cosine(w, ez)

    # make all of the above into a 3 X 3 matrix and return it
    c_matrix = np.asarray(
        [
            [u_alpha, u_beta, u_gamma],
            [v_alpha, v_beta, v_gamma],
            [w_alpha, w_beta, w_gamma],
        ]
    )

    return c_matrix
