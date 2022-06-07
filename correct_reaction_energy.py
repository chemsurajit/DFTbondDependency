import sys
import logging
import argparse
import xyz2mol


eV2au = 0.03674930495120813 # Hartree
eV2kcal = 23.060541945329334 # kcal/mole


bonds_correction_pbe = {'C_d_C': -0.04720942012121344,
                        'C_t_C': -0.21051127696657002,
                        'C_C_A': 0.01350706283940116,
                        'C_s_O': 0.058913661477387906,
                        'C_d_O': 0.09478783968914134,
                        'C_O_A': 0.058601727073329286,
                        'C_s_N': -0.01700493870792873,
                        'C_d_N': -0.04393126692685779,
                        'C_t_N': -0.3335535240446993,
                        'C_N_A': -0.02169970912158023,
                        'O_s_N': 0.17617543138309522,
                        'O_d_N': 0.2982693293763882,
                        'O_N_A': 0.1336410237248417,
                        'N_s_N': 0.19291401513188064,
                        'N_d_N': 0.007537011973101717,
                        'N_N_A': 0.011176795612685095}


bonds_correction_b3lypd = {'C_d_C': 0.11002126797437362,
                           'C_t_C': 0.12708094409389806,
                           'C_C_A': 0.06514956470815053,
                           'C_s_O': 0.07810842075166526,
                           'C_d_O': 0.2553188263219433,
                           'C_O_A': 0.08926953833568295,
                           'C_s_N': -0.009121212284876754,
                           'C_d_N': 0.0999968172590848,
                           'C_t_N': 0.11670813250640646,
                           'C_N_A': 0.03254948736725622,
                           'O_s_N': 0.08517758567921756,
                           'O_d_N': 0.23783366526756422,
                           'O_N_A': 0.11788819117347499,
                           'N_s_N': 0.06530770276077955,
                           'N_d_N': 0.09504374302205423,
                           'N_N_A': 0.014127992860272899}


bonds_correction_m062x = {'C_d_C': -0.04875042692605581,
                          'C_t_C': -0.08714341858770626,
                          'C_C_A': -0.004005408806382787,
                          'C_s_O': 0.04143671641887671,
                          'C_d_O': -0.03503331595948148,
                          'C_O_A': 0.02564221120767289,
                          'C_s_N': -0.008323388498284717,
                          'C_d_N': -0.07731454017306497,
                          'C_t_N': -0.12312842822764032,
                          'C_N_A': -0.015247188823177648,
                          'O_s_N': -0.06304967509565461,
                          'O_d_N': -0.21819014074960713,
                          'O_N_A': -0.08968208415965984,
                          'N_s_N': -0.028201960939824758,
                          'N_d_N': -0.0700270847383631,
                          'N_N_A': -0.07342594474772787}


bonds_names_list = [k for k, v in bonds_correction_pbe.items()]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="To calculate the energy correction using the LR results."
    )
    parser.add_argument(
        "-eunit", "--eunit",
        type=str,
        default="eV",
        choices=["eV", "a.u.", "kcal/mol", "kJ/mol"],
        help="Energy unit to be included."
    )
    parser.add_argument(
        "-r", "--r",
        type=str,
        required=True,
        nargs="+",
        default=None,
        help="Reactant log files in ADF format"
    )
    parser.add_argument(
        "-p", "--p",
        type=str,
        required=True,
        nargs="+",
        default=None,
        help="Product log files in ADF format."
    )
    parser.add_argument(
        "-dft_functional", "--dft_functional",
        type=str,
        default=None,
        required=True,
        choices=["PBE", "B3LYP-D", "M06-2X"],
        help="Name of the DFT functional used to compute the energies"
    )
    parser.add_argument(
        "-postscf", "--postscf",
        required=False,
        action="store_true",
        help="Switch that indicates whether the energies were computed with postSCF or not."
    )
    return parser.parse_args()


def get_coords_frm_adflog(logfile):
    line_nos = []
    atoms = []
    coords = []
    with open(logfile) as flog:
        for line_no, line in enumerate(flog, 1):
            if line.strip() == "Atoms":
                line_nos.append(line_no)

    with open(logfile) as flog:
        last_coords_to_end = flog.readlines()[(line_nos[-1]+1):]
        for line in last_coords_to_end:
            splitline = line.split()
            if len(splitline) == 5:
                atoms.append(splitline[1])
                icoord = [
                    float(splitline[2]),
                    float(splitline[3]),
                    float(splitline[4])
                ]
                coords.append(icoord)
            else:
                break
    assert len(atoms) == len(coords)
    return atoms, coords


def get_energy(logfile, eunit, dft_functional, post_scf):
    len_content = ""
    energy = 0.0
    with open(logfile) as fp:
        if post_scf:
            en_dict = {}
            for lines in fp.readlines():
                if "FR:" in lines:
                    functional = lines[4:19].strip().upper()
                    if eunit == "a.u.":
                        fid = 0
                    elif eunit == "eV":
                        fid = 1
                    elif eunit == "kcal/mol":
                        fid = 2
                    elif eunit == "kJ/mol":
                        fid = 3
                    en = float(lines.split("=")[1].split()[fid])
                    en_dict[functional] = en
            if dft_functional in en_dict:
                energy = en_dict[dft_functional]
            else:
                logging.error("Energy for postSCF DFT (%s) not found in %s" % (dft_functional, logging))
        else:
            for lines in fp.readlines():
                if "Total Bonding Energy:" in lines:
                    len_content = lines
                if len_content:
                    splitted_line = len_content.split()
                    if eunit == "a.u.":
                        energy = float(splitted_line[3])
                    elif eunit == "eV":
                        energy = float(splitted_line[4])
                    elif eunit == "kcal/mol":
                        energy = float(splitted_line[5])
                    elif eunit == "kJ/mol":
                        energy = float(splitted_line[6])
                    else:
                        logging.error("Unit not recognized: %s" % eunit)
                else:
                    logging.error("Log file: %s doesn't contain energy." % logfile)
    return energy


def get_mol_data_from_adf_logs(logfiles, post_scf, dft_functional, eunit=None):
    if eunit is None:
        eunit = "eV"
    bonds_count_dict = {bonds: 0 for bonds in bonds_names_list}
    dft_energy = 0.0
    for logfile in logfiles:
        atoms, coords = get_coords_frm_adflog(logfile)
        atoms = [xyz2mol.int_atom(atom) for atom in atoms]
        dft_energy += get_energy(logfile, eunit, dft_functional, post_scf)
        mols = None
        logging.debug("atoms: %s" % atoms)
        logging.debug("coords: %s" % coords)
        try:
            mols = xyz2mol.xyz2mol(atoms, coords, charge=0.0,
                               use_graph=True, allow_charged_fragments=False,
                               embed_chiral=True, use_huckel=False)
        except:
            logging.warning("Failed to convert mol object for file %s" % logfile)

        if mols is not None:
            for bond in mols[0].GetBonds():
                # RDkit can consider a C=O bond (eg) as either C=O or O=C.
                # So, to compare the strings, bonds from both the front side
                # i.e., from C side and back side, i.e., O side is considered.
                bname_front = None
                bname_back = None
                bond_type = bond.GetBondType()
                atom1 = mols[0].GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()
                atom2 = mols[0].GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()
                if str(bond_type) == "SINGLE":
                    bname_front = atom1 + "_" + "s_" + atom2
                    bname_back = atom2 + "_" + "s_" + atom1
                if str(bond_type) == "DOUBLE":
                    bname_front = atom1 + "_" + "d_" + atom2
                    bname_back = atom2 + "_" + "d_" + atom1
                if str(bond_type) == "TRIPLE":
                    bname_front = atom1 + "_" + "t_" + atom2
                    bname_back = atom2 + "_" + "t_" + atom1
                if str(bond_type) == "AROMATIC":
                    bname_front = atom1 + "_" + atom2 + "_" + "A"
                    bname_back = atom2 + "_" + atom1 + "_" + "A"
                if str(bond_type) not in ['SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC']:
                    logging.warning("Unknown bond type: ", bond_type, " Will not provide correction")
                if bname_front in bonds_count_dict:
                    bonds_count_dict[bname_front] += 1
                elif bname_back in bonds_count_dict:
                    bonds_count_dict[bname_back] += 1
                else:
                    logging.warning(
                        "No correction term for: %s or %s" % (bname_front, bname_back)
                    )
    return dft_energy, bonds_count_dict


def get_bond_correction(reactant_bond_dict, pdt_bond_dict, lr_coeff_dict, eunit):
    total_correction = 0.0
    for k, v in reactant_bond_dict.items():
        if eunit == "kcal/mol":
            v = v*eV2kcal
        if eunit == "a.u.":
            v = v*eV2au
        total_correction -= v*lr_coeff_dict[k]
    for k, v in pdt_bond_dict.items():
        if eunit == "kcal/mol":
            v = v*eV2kcal
        if eunit == "a.u.":
            v = v*eV2au
        total_correction += v*lr_coeff_dict[k]
    return total_correction


def main():
    args = parse_arguments()
    reactant_log_files = args.r
    product_log_files = args.p
    eunit = args.eunit
    post_scf = args.postscf
    dft_functional = args.dft_functional
    reactant_energy, reactant_bond_dict = get_mol_data_from_adf_logs(
        reactant_log_files, post_scf, dft_functional, eunit=eunit
    )
    pdt_energy, pdt_bond_dict = get_mol_data_from_adf_logs(
        product_log_files, post_scf, dft_functional, eunit=eunit
    )
    lr_coeff_dict = {}
    if dft_functional == "PBE":
        lr_coeff_dict = bonds_correction_pbe
    elif dft_functional == "B3LYP-D":
        lr_coeff_dict = bonds_correction_b3lypd
    elif dft_functional == "M06-2X":
        lr_coeff_dict = bonds_correction_m062x
    reaction_energy = pdt_energy - reactant_energy
    total_bond_correction = get_bond_correction(
        reactant_bond_dict, pdt_bond_dict, lr_coeff_dict, eunit
    )
    print("%s energy: %f %s" % (dft_functional, reaction_energy, eunit))
    print("%s corrected energy: %f %s" % (dft_functional, (reaction_energy + total_bond_correction), eunit))
    print("Correction to %s: %f %s" % (dft_functional, total_bond_correction, eunit))
    return


if __name__ == "__main__":
    main()
