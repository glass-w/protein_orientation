def calc_lipid_contacts(universe, cut_off):

        # reset the centre of mass list for each residue
        com_array = []

        # Cycle through every residue and get it's centre of mass based on it's side chain atoms (from dict)
        for resid in resid_list:
            sel = universe.select_atoms("resid " + str(resid))

            com_array.append(sel.select_atoms("name " + str(term_atom_dict[sel.resnames[0]])).atoms.center_of_mass())

        # convert centre of mass array (a list) for next step, define the head group of the POPC groups.
        com_array = np.array(com_array, dtype="float32")

        popc_head = universe.select_atoms("resname POPC and name P")

        if cut_off:
            # Calc distances of CoMs of residues and POPC heads within a cut off, returns a boolean array
            d = dist(com_array, popc_head.atoms.positions) <= cut_off

        else:
            # Calc distances of CoMs of residues and POPC heads
            d = dist(com_array, popc_head.atoms.positions)

        return d
