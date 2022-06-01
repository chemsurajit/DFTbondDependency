import os
import argparse
import xyz2mol
from collections import Counter
from scm import plams
from pathlib import Path


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
    parser = argparse.ArgumentParser(description="To calculate the energy correction using the LR results.")
    parser.add_argument('--eunit',
                        type=str,
                        default='eV',
                        choices=['eV', 'a.u.', 'kcal/mol'],
                        help="Energy unit to be included.")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--react_xyz_files',
                        type=str,
                        nargs='+',
                        default=None,
                        help="Reactant xyz files")
    requiredNamed.add_argument('--react_out_file',
                        type=str,
                        nargs='+',
                        default=None,
                        help="Reactant log/out file.")
    requiredNamed.add_argument('--pdt_xyz_files',
                        type=str,
                        nargs='+',
                        default=None,
                        help="Product xyz files.")
    requiredNamed.add_argument('--pdt_out_file',
                        type=str,
                        nargs='+',
                        default=None,
                        help="Product log/out file.")
    requiredNamed.add_argument('--dft_functional',
                        type=str,
                        default=None,
                        choices=["PBE", "B3LYP", "M062X"],
                        help="Name of the DFT functional used to compute the energies")
    return parser.parse_args()


def get_bonds_count_from_xyz(xyzfile=None):
    bonds_count_dict = {bonds: 0 for bonds in bonds_names_list}
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(xyzfile)
    mols = None
    try:
        mols = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge=0,
                               use_graph=True, allow_charged_fragments=False,
                               embed_chiral=True, use_huckel=False)
    except:
        print("No mol object from the xyz file: %s " % xyzfile)
        raise Warning("xyzfile not found")
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
            #print("debug: ", bond_type, atom1, atom2)
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
                print("Unknown bond type: ", bond_type, " Will not provide correction")
            if bname_front in bonds_count_dict:
                bonds_count_dict[bname_front] += 1
            elif bname_back in bonds_count_dict:
                bonds_count_dict[bname_back] += 1
            else:
                print("Bonds either: %s or: %s does not exists in the bonds list" % (bname_front, bname_back))
    else:
        raise Warning("RDKit mol object not found. Conversion problem")
    return bonds_count_dict


def get_dft_energy(outfile, unit='eV'):
    """This function returns the energy in hartree and ev unit. Implementation only for SCM.
       outdir is the directory where scm job ran."""
    energy = None
    with open(outfile) as fp:
        for lines in fp.readlines():
            if "Bond Energy" in lines:
                line_content = lines.split()
                energy = float(line_content[-2])
                enunit = line_content[-1]
                if enunit == unit:
                    break
    return energy


def main():
    args = parse_arguments()
    print("List of bonds to be corrected: ", bonds_names_list)
    # first get the bond list and the energies for each of the reactants and products.
    dft_energy = 0.0
    reaction_bond_change_dict = {bond: 0 for bond in bonds_names_list}
    nreactant = len(args.react_xyz_files)
    eunit = args.eunit
    for i in range(nreactant):
        react_xyz = args.react_xyz_files[i]
        react_out_file = args.react_out_file[i]
        for k, v in get_bonds_count_from_xyz(react_xyz).items():
            reaction_bond_change_dict[k] -= v
        # for the reactant side, the energy will be subtracted.
        dft_energy -= get_dft_energy(react_out_file, unit=eunit)
        print("reactant: ", reaction_bond_change_dict)
    npdt = len(args.pdt_xyz_files)
    for i in range(npdt):
        pdt_xyz = args.pdt_xyz_files[i]
        pdt_out_file = args.pdt_out_file[i]
        for k,v in get_bonds_count_from_xyz(pdt_xyz).items():
            reaction_bond_change_dict[k] += v
        print("product: ", reaction_bond_change_dict)
        # for the product side, the energy will be added.
        dft_energy += get_dft_energy(pdt_out_file, unit=eunit)
    # correction to be done from the coefficients
    print("DFT reaction energy: %f %s" %(dft_energy, eunit))
    print("bond changes: ", reaction_bond_change_dict)
    if args.dft_functional == "PBE":
        correction_dict = bonds_correction_pbe
    elif args.dft_functional == "B3LYP":
        correction_dict = bonds_correction_b3lypd
    elif args.dft_functional == "M062X":
        correction_dict = bonds_correction_m062x
    for k, v in reaction_bond_change_dict.items():
        if eunit == "kcal/mol":
            v = v*eV2kcal
        if eunit == "a.u.":
            v = v*eV2au
        dft_energy += v * correction_dict[k]
    print("Error corrected reaction energy: %f %s" %(dft_energy, eunit))
    return


if __name__ == "__main__":
    main()
