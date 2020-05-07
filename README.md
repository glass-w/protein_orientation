[![Build Status](https://travis-ci.org/WG150/protein_orientation.svg?branch=refactor)](https://travis-ci.org/WG150/protein_orientation)
[![codecov](https://codecov.io/gh/WG150/protein_orientation/branch/refactor/graph/badge.svg)](https://codecov.io/gh/WG150/protein_orientation)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

protein_orientation.py is a script to analyse the pitch of a protein region over the course of a molecular dynamics simulation. 

## Note

Although this script has been made available it has only been tested for biomolecular systems specific to one research project. As a result the code **is not guaranteed to work with all systems and / or provide the correct result**. Please use this code at your own risk. Please raise any issues you encounter and contribute if you can.  

## Usage
  
```bash
python protein_orientation.py -h
```
```text
usage: protein_orientation.py [-h] -c GRO_FILE_LIST -f XTC_FILE_LIST -com_sel
                              COM_SELECTION [-method METHOD] -n
                              NUM_OF_PROTEINS [-skip SKIP]
                              [-vtraj VECTOR_TRAJ]
                              [-res_vector_sel RES_VECTOR_SEL]
                              [-stride STRIDE_FILE] [-pa_only PA_SINGLE]
                              [-nprocs NPROCS] [-ref_option REF_OPTION]
                              [-sec_struc_choice SEC_STRUC_CHOICE]

Calculates the orientation of a user defined region of a protein

optional arguments:
  -h, --help            show this help message and exit
  
  -c GRO_FILE_LIST      The list of coordinate files [.gro], this takes the
                        form of a text file with each file location starting
                        on a new line.
  
  -f XTC_FILE_LIST      The list of corrected trajectory files: pbc artifacts
                        removed, no jumping across PBC. This takes the form of
                        a text file with each file location starting on a new
                        line.
  
  -com_sel COM_SELECTION
                        The range of resids to use for centre of mass
                        calculation, in the form of A:B, where A and B are
                        integers.
  
  -method METHOD        The vectors can be calculated by 1) a set of user
                        defined vectors based on the centre of mass of the
                        main selection and the alpha carbon (CA) of a
                        specified residue OR 2) the method can be used in
                        combination with (1) and use the principal axes of
                        inertia. In either (1) or (2) the user must define a
                        set of vectors that roughly correspond to the
                        principal axes - this ensures that when calculated
                        they always point in the direction specified by the
                        users vectors. Options: user or user_pa. Default =
                        user_pa
  
  -n NUM_OF_PROTEINS    Number of protein copies in the system, default 1.
  
  -skip SKIP            The number of frames to skip, default 1.
  
  -vtraj VECTOR_TRAJ    Set to True if you want a trajectory of the vectors,
                        default False.
  
  -res_vector_sel RES_VECTOR_SEL
                        The resids of the residues to use for the roll, pitch,
                        and yaw calculation respectively: in the form A, B, C.
  
  -stride STRIDE_FILE   The name of the stride file to read, a .txt file. This
                        will be used in combination with the -com_sel
                        selection to only choose those residues involved in
                        secondary structure. If using the 'user_pa' method
                        (see -method) this option must be supplied.
  
  -pa_only PA_SINGLE    If set to True a principal component calculation will
                        be carried out and written to a .pdb file, this is to
                        help in selecting the appropriate residues for a run.
                        Default False
  
  -nprocs NPROCS        Number of processes to use, default=1.
  
  -ref_option REF_OPTION
                        Choice of what basis of vectors to use as a reference,
                        from which the Euler angles will be calcualted.
                        Permitted chocies are: 'first_frame', angles will be
                        calculated in reference to the PAs calculated in the
                        first frame. 'user', angles will be calculated in
                        reference to a user defined set of vectors. 'standard'
                        where the standard is x, y, z = [1,0,0], [0,1,0],
                        [0,0,1]. default = 'standard'.

  -ref_basis REF_BASIS  To be used in combination with 'user' if used (see
                        -ref_option). The basis vectors to be used as a
                        reference, if not passed the default will be used (see
                        -ref_option). This should be a .txt file with the x,
                        y, z coordinates on each line.  

  -sec_struc_choice SEC_STRUC_CHOICE
                        A file containing the choice of secondary structure to
                        use in the calculation of the centre of mass. If using
                        the 'user_pa' method (see -method) this option must be
                        supplied. Valid choices include: 'strand', '310helix',
                        or 'alphahelix'. In the file these must be comma
                        separated and have no whitespace between them. e.g.
                        strand,310helix
```

## Example usage

```text
python orientation/protein_orientation.py -c data/gro_file_list.txt -f data/xtc_file_list.txt -com_sel 1:123 -n 1 -method user_pa -res_vector_sel 102,69,5 -stride data/beta3_stride_file.txt -nprocs 2 -skip 10 -ref_option standard -sec_struc_choice data/sec_struc.txt
```

Tests
------
If you would like to run written tests and / or calculate code coverage use the following commands:
```text
python -m pytest ./
python -m pytest --cov=./
```

Thanks
------
Thanks go to Rocco Meli, Aphroditi Zaki, and Irfan Alibay for mathematical discussions and support in code development.

Notes
-----
* This code was used to calculate protein orientation in [this publication](https://www.frontiersin.org/articles/10.3389/fmolb.2020.00040/full).
