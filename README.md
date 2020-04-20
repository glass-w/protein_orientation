# Protein Orientation

This is the inital version of protien_orientation - currently undergoing refactoring (see branch)

## Usage
  
```bash
python protein_orientation.py -h
```

```text
usage: protein_orientation.py [-h] [-c GRO_FILE_LIST] [-f XTC_FILE_LIST]
                              [-com_sel COM_SELECTION] [-n NUM_OF_PROTEINS]
                              [-skip SKIP] [-vtraj VECTOR_TRAJ]
                              [-method METHOD]
                              [-res_vector_sel RES_VECTOR_SEL]
                              [-stride STRIDE_FILE] [-pa_only PA_SINGLE]
                              [-nprocs NPROCS] [-ref_option REF_OPTION]
                              [-ref_basis REF_BASIS]
                              [-sec_struc_choice SEC_STRUC_CHOICE]

optional arguments:
  -h, --help            show this help message and exit
  -c GRO_FILE_LIST      The list of coordinate files [.gro], this takes the
                        form of a text file
  -f XTC_FILE_LIST      A list of corrected trajectory files: pbc artifacts
                        removed, no jumping across PBC, thisis a text file.
  -com_sel COM_SELECTION
                        The range of resids to use for centre of mass
                        calculation, in the form of A:B where A and B are
                        integers
  -n NUM_OF_PROTEINS    The number of protein copies in the system
  -skip SKIP            The number of frames to skip
  -vtraj VECTOR_TRAJ    Set to True if you want a trajectory of the vectors
  -method METHOD        The vectors can be calculated by 1) a set of user
                        defined vectors based on the centre of mass of the
                        main selection and the CA of a specified residue OR 2)
                        the method can be used in combination with (1) and use
                        the principal axes of inertia. In either (1) or (2)
                        the user must define a set of vectors that roughly
                        correspond to the principal axes - this ensures that
                        when calculated they always point in the direction
                        specified by the user's vectors Options: 'user'
                        'user_pa'
  -res_vector_sel RES_VECTOR_SEL
                        The resids of the residues to use for the roll, pitch,
                        and yaw calculation in the form A, B, C.
  -stride STRIDE_FILE   the name of the stride file to read, a .txt file. This
                        will be used in combination with the -com_sel
                        selection to only choose those residues involved in
                        secondary structure.
  -pa_only PA_SINGLE    If set to True a principal component calculation will
                        be carried out and written to a .pdb file, this is to
                        help in selecting the appropriate residues for a run.
  -nprocs NPROCS        No. or processes to use, usually set to the number of
                        repeats you have. Max = max No. of CPUs available to
                        you
  -ref_option REF_OPTION
                        Choice of what basis of vectors to use as a reference,
                        from which the Euler angles will be calcualted.
                        Permitted chocies are: 'first_frame', 'user' or
                        'standard' where standard is x, y, z = [1,0,0],
                        [0,1,0], [0,0,1].
  -ref_basis REF_BASIS  The basis vectors to be used as a reference, if not
                        passed the default will be used (see -ref_option).
                        This should be a .txt file with the x, y, z
                        coordinates one each line e.g. 1,0,0 0,1,0 0,0,1
  -sec_struc_choice SEC_STRUC_CHOICE
                        A file containing the choice of secondary srtucture to
                        use in the calculation of the centre of mass. Valid
                        choices include: strand, 310helix, or alphahelix. In
                        the file these must be comma separated and have no
                        whitespace between them. e.g. strand,310helix

```

