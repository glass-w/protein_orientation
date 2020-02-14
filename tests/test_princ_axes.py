import MDAnalysis as mda

from orientation import calc_angles

u = mda.Universe('data/b3_syst_protein_only.gro')

def check_pas():

	# select alpha carbons only
    CA = universe.select_atoms(selection)

    # calculate moment of inertia, and sort eigenvectors
    I = CA.moment_of_inertia()

    # UT = CA.principal_axes()
    # U = UT.T

    values, evecs = np.linalg.eigh(I)
    indices = np.argsort(values)
    U = evecs[:, indices]

    Lambda = U.T.dot(I.dot(U))
    
    print(Lambda)
    print(np.allclose(Lambda - np.diag(np.diagonal(Lambda)), 0))

    print("")
    print(U) 
