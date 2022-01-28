import argparse
import os
import sys
import fnmatch
import xyz2mol
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import pandas as pd
import csv


def get_files(directory, match=""):
    """This function will return a list of xyz files as list. The list elements
    will be in the form of path object. The expected file name is: dsgdb9nsd_n.xyz
    where n is an integer number.
    args:
        directory - Path from which xyz files needs to be found.
        match - a matching string for the files. Else, all the files will be returned.
    returns:
        xyz_files - List of xyz files with names dsgdb9nsd_*.xyz after being sorted."""
    allfiles = []
    for files in os.listdir(directory):
        # Take the xyz files which has 'dsgdb9nsd_' in their name
        if fnmatch.fnmatch(files, match):
            allfiles.append(os.path.abspath(os.path.join(directory, files)))
    # Check if the list xyz_files is empty. If it is empty, then program exits.
    if not allfiles:
        print("No files with match %s found in directory: %s " % (match, directory))
        print("Program exit now.")
        sys.exit()
    return sorted(allfiles)


def get_dft_energies(logfile):
    """This function takes a ADF logfile and return the dft energies as dictionary.
    args:
        logfile - ADF output file
    returns:
        dft_energies - dictionary of DFT energies. Keys as functional, value as energy in eV.
    """
    dft_energies = {}
    with open(logfile, 'r') as flog:
        for line in flog:
            if "FR:" in line:
                FR_cont = True
                func = line[4:19].strip().upper()
                en_ev = float(line.split("=")[1].split()[1])
                dft_energies[func] = en_ev
    if not FR_cont:
        print ("No data found for DFT functional")
        print ("The reason could be that the ADF log file prints the energy with different format.")
    return dft_energies


def mol_to_bonds_list(molobj, filename):
    """This function takes a RDKit mol object as input with xyzfile and makes
    a list of bonds. The list of bonds will always be confined in the variable
    defined as bonds_list. If some bonds not present, it will be assigned to zero.
    args:
        molobj - RDKit format mol object
        filename - The xyz file name
    returns:
        bonds_list - list of bonds after counting each of the bond types.
    """
    bonds_list = {'C_s_C': 0, 'C_d_C': 0, 'C_t_C': 0, 'C_C_A': 0,
                  'C_s_H': 0, 'C_s_O': 0, 'C_d_O': 0, 'C_t_O': 0,
                  'C_O_A': 0, 'C_s_N': 0, 'C_d_N': 0, 'C_t_N': 0,
                  'C_N_A': 0, 'C_s_F': 0, 'O_s_O': 0, 'O_d_O': 0,
                  'O_O_A': 0, 'O_s_H': 0, 'O_s_N': 0, 'O_d_N': 0,
                  'O_t_N': 0, 'O_N_A': 0, 'O_s_F': 0, 'N_s_N': 0,
                  'N_d_N': 0, 'N_t_N': 0, 'N_N_A': 0, 'N_s_H': 0,
                  'N_s_F': 0, 'F_s_H': 0
                  }

    xyzfile = os.path.basename(filename)
    if molobj is not None:
        for bond in molobj.GetBonds():
            btype = bond.GetBondType()
            a1 = molobj.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()
            a2 = molobj.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()
            #
            # C-C bonds
            #
            if a1+a2 == "CC":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_C'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_C'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_C'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_C_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), a1+"-"+a2)
                    print("Exiting for.", xyzfile)
                    raise ValueError
            #
            # C-H bonds
            #
            elif a1+a2 == "CH" or a2+a1 == "CH":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_H'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # C-O bonds
            #
            elif a1+a2 == "CO" or a2+a1 == "CO":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_O'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_O'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_O'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_O_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # C-N bonds
            #
            elif a1+a2 == "CN" or a2+a1 == "CN":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_N_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # C-F bonds
            #
            elif a1+a2 == "CF" or a2+a1 == "CF":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_F'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # O-O bonds
            #
            elif a1+a2 == "OO":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_O'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['O_d_O'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['O_O_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # O-H bond
            #
            elif a1+a2 == "OH" or a2+a1 == "OH":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_H'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # O-N bonds
            #
            elif a1+a2 == "ON" or a2+a1 == "ON":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['O_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['O_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['O_N_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # O-F bond
            #
            elif a1+a2 == "OF" or a2+a1 == "OF":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_F'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # N-N bonds
            #
            elif a1+a2 == "NN":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['N_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['N_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['N_N_A'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # N-H bond
            #
            elif a1+a2 == "NH" or a2+a1 == "NH":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_H'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # N-F bonds
            #
            elif a1+a2 == "NF" or a2+a1 == "NF":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_F'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
            #
            # F-H bond
            #
            elif a1+a2 == "FH" or a2+a1 == "FH":
                if str(btype) == "SINGLE":
                    bonds_list['F_s_H'] += 1
                else:
                    print("Don't know bond type: ", str(btype), " ", a1+"-"+a2)
                    print("Exiting for ", xyzfile)
                    raise ValueError
    return bonds_list


def get_properties_combined(index, smiles, chem_formula, bonds_list):
    """This function takes an index, smiles string and a list of bonds and return
    a dictionary with the proper key names. This key names will be used homogeneously
    in all the programs that follows.
    args:
        index - index of the molecule in the QM9_GMP2 dataset
        smiles - smiles string as calculated using the xyz2mol program
        chem_formula - Formula as type string.
        bonds_list - A dictionary with all the bond names as keys and their numbers
                     in the molecule as values.
    returns:
        dictionary with keys as below.
    """
    all_properties = {"index": index, "smiles": smiles, "chemformula":chem_formula}
    return {**all_properties, **bonds_list}


def get_smiles_from_xyz(ixyzfile):
    """Function that takes xyzfile path and return the smiles"""
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(ixyzfile)
    return_code = True
    mols = None
    try:
        mols = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge=0,
                               use_graph=True, allow_charged_fragments=False,
                               embed_chiral=True, use_huckel=False)
    except:
        mols = None
        print("No mol object from the xyz file: ", ixyzfile)
    return mols


def update_failed_indices(outputfile, indices):
    if not indices:
        print("Number of failed indices = ", len(indices))
        print("No file will be created.")
    else:
        with open(outputfile, 'w') as fp:
            for index in indices:
                fp.write("%s\n" % str(index))
    return


def update_pd_df(inpdf, index=None, smiles=None, chemformula=None, bonds=None, g4mp2_energy=None, energies=None):
    """Function to update to dataframe
    """
    row_dict = {"index":index, "smiles":smiles, "chemformula":chemformula, "G4MP2":g4mp2_energy}
    row_dict.update(energies)
    row_dict.update(bonds)
    df = pd.DataFrame(row_dict, index=[0])
    newdf = pd.concat([inpdf, df])
    return newdf


def main():
    #
    # parsing arguments
    #
    parser = argparse.ArgumentParser("File to create bond list from xyz file.")
    parser.add_argument('-xd', '--xyzd',
                        help="Location of the xyz directory. Default is current directory.",
                        required=True)
    parser.add_argument('-ld', '--dir_logfiles',
                        help="Location of the directory from where the logfiles will be read.",
                        required=True)
    parser.add_argument('-g', '--g4mp2_csvfile',
                        help="Location of the csv file containing G4MP2 energies of the QM9 dataset.", required=True)
    parser.add_argument('-o', '--output', help="Name of the output file. Should end with csv",
                        required=False, default="qm9_bonds_energies.csv")
    parser.add_argument('-f', '--failed_file', help="File in which to update the indices where processing failed",
                        required=False, default='failed_indices.dat')


    args = parser.parse_args()
    xyz_direcoty = os.path.abspath(args.xyzd)
    log_dir = args.dir_logfiles
    qm9_g4pm2_csv = args.g4mp2_csvfile
    output_csv = args.output
    failed_output = args.failed_file
    #
    # End of parsing input arguments
    #
    failed_mols_indices = []
    print("Directory containing xyz files: ", xyz_direcoty)
    xyzfiles = get_files(xyz_direcoty, match='dsgdb9nsd_*.xyz')
    logfiles = get_files(log_dir, match='*_xyz.out')
    print("Number of xyz files: ", len(xyzfiles))
    df = pd.DataFrame()
    g4mp2_pd = pd.read_csv(qm9_g4pm2_csv)
    for i in range(len(xyzfiles)):
        ixyzfile = xyzfiles[i]
        ilogfile = logfiles[i]
        logindex = os.path.basename(ilogfile).split("_")[0]
        xyzindex = os.path.basename(ixyzfile).split("_")[1].split(".")[0]
        if logindex != xyzindex:
            print("The index in logfile and xyzfiles are not same.")
            print("Xyzfile: ", ixyzfile)
            print("logfile: ", ilogfile)
            break
        print(logindex, ilogfile, ixyzfile)
        dft_energies = get_dft_energies(os.path.abspath(ilogfile))
        if not dft_energies:
            print("No dft energies found for index: ", logindex)
            failed_mols_indices.append(logindex)
            continue
            # Get the smile string using the xyz2mol module.
        mols = get_smiles_from_xyz(ixyzfile)
        if mols is None:
            print("Failed to convert xyz to smiles string: ", ixyzfile)
            failed_mols_indices.append(xyzindex)
            continue
        bonds_list = mol_to_bonds_list(mols[0], ixyzfile)
        smiles = Chem.MolToSmiles(mols[0], isomericSmiles=True)
        chemformula = CalcMolFormula(mols[0])
        g4mp2_energy = g4mp2_pd.loc[g4mp2_pd['index'] == int(logindex), 'G4MP2'].tolist()[0]
        #
        # Write everything to csv file
        #
        df = update_pd_df(df, index=logindex, smiles=smiles, chemformula=chemformula, bonds=bonds_list, g4mp2_energy=g4mp2_energy, energies=dft_energies)
        if (i % 10000) == 0:
            print("Number of molecules converted: ", i)
    #write to csv file.
    df.to_csv(output_csv, sep=",", index=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
    update_failed_indices(failed_output, failed_mols_indices)
    return


if __name__ == "__main__":
    main()
