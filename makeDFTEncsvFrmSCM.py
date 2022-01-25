import os, sys
import fnmatch
import csv
import pandas as pd
import argparse
import xyz2mol
from rdkit import Chem


def get_all_logfiles(log_dir):
    """This function will return a list of all the ADF logfiles. The name of
    the logfiles should be like: n_xyz.out where n is integer.
    args:
        log_dir - The location where the logfiles are kept
    returns:
        logfiles - a sorted list of all the log file names (with path)"""
    logfiles = []
    for files in os.listdir(log_dir):
        if fnmatch.fnmatch(files, '*_xyz.out'):
            logfiles.append(os.path.abspath(os.path.join(log_dir, files)))
    # Check if the list logfiles is empty:
    if not logfiles:
        print("No files with name *_.out found in directory: ", log_dir)
        print("Program exit now.")
        sys.exit()
    return sorted(logfiles)


def get_xyz_files(xyz_dir):
    """This function will return a list of xyz files as list. The list elements
    will be in the form of path object. The expected file name is: dsgdb9nsd_n.xyz
    where n is an integer number.
    args:
        directory - Path from which xyz files needs to be found.
    returns:
        xyz_files - List of xyz files with names dsgdb9nsd_*.xyz after being sorted."""
    xyz_files = []
    for files in os.listdir(xyz_dir):
        # Take the xyz files which has 'dsgdb9nsd_' in their name
        if fnmatch.fnmatch(files, 'dsgdb9nsd_*.xyz'):
            xyz_files.append(os.path.abspath(os.path.join(xyz_dir, files)))
    # Check if the list xyz_files is empty. If it is empty, then program exits.
    if not xyz_files:
        print("No files with extension .xyz found in directory: ", xyz_dir)
        print("Program exit now.")
        sys.exit()
    return sorted(xyz_files)


def save_to_output(outputfile, write_mode, key_value_pairs):
    """Function to write to csv file"""
    df = pd.DataFrame(key_value_pairs, index=[0])
    if write_mode.lower() == "w":
        df.to_csv(outputfile, mode=write_mode, sep=",", index=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
    elif write_mode.lower() == "a":
        df.to_csv(outputfile, mode=write_mode, sep=",", index=False, header=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
    else:
        print("CSV file operating mode", write_mode, " not understood")
    return

def get_dft_energies(logfile):
    """This function takes a ADF logfile and return the dft energies as dictionary.
    args
        logfile - ADF output file
    returns
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


def get_smiles_from_xyz(ixyzfile):
    """Function that takes xyzfile path and return the smiles"""
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(ixyzfile)
    return_code = True
    try:
        mols = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge=0,
                           use_graph=True, allow_charged_fragments=False,
                           embed_chiral=True, use_huckel=False)
    except:
        return_code = False
    if return_code:
        if len(mols) == 0:
            print("No mol object from file: ", files)
        smiles = Chem.MolToSmiles(mols[0], isomericSmiles=True)
    else:
        smiles = ""
    return smiles


def update_failed_indices(outputfile, indices):
    with open(outputfile, 'w') as fp:
        for index in indices:
            fp.write("%s\n" % str(index))
    return


def main():
    parser = argparse.ArgumentParser("Program to create csv files from the logfiles of SCM output.")
    parser.add_argument('-d', '--dir_logfiles',
                        help="Location of the directory from where the logfiles will be read.",
                        required=False, default='./')
    parser.add_argument('-xd', '--xyzd', help="Directory where the xyz files are present.",
                        required=True)
    parser.add_argument('-o', '--output_csv', help="Name of the output csv file for the energies.",
                        required=False, default='qm9_dft_en.csv')
    parser.add_argument('-f', '--failed_file', help="File in which to update the indices where processing failed",
                        required=False, default='failed_dft.dat')


    args = parser.parse_args()
    xyz_dir = args.xyzd
    output_csv = args.output_csv
    log_dir = args.dir_logfiles
    failed_output = args.failed_file


    all_key_vals = {}
    print("XYZ dir: ", xyz_dir)
    print("output csv: ", output_csv)
    print("log directory: ", log_dir)
    logfiles = get_all_logfiles(log_dir)
    print("Number of logfiles: ", len(logfiles))
    xyzfiles = get_xyz_files(xyz_dir)
    print("Number of xyzfiles: ", len(xyzfiles))
    if len(logfiles) != len(xyzfiles):
        print("The number of logfiles and xyzfiles are different. Please check the directories.")
        print("The program will exit now.")
        sys.exit()


    failed_file_indices = []


    for i in range(len(logfiles)):
        if i == 0:
            write_mode = 'w'
        else:
            write_mode = 'a'
        ilogfile = logfiles[i]
        ixyzfile = xyzfiles[i]
        # Then, check if the indices of the logfiles and xyzfiles are same
        # This is a kind of sanity check
        logindex = os.path.basename(ilogfile).split("_")[0]
        xyzindex = os.path.basename(ixyzfile).split("_")[1].split(".")[0]
        if logindex != xyzindex:
          print("The index in logfile and xyzfiles are not same.")
          print("Xyzfile: ", ixyzfile)
          print("logfile: ", ilogfile)
          break
        print(logindex, ilogfile, xyzfiles[i])
        # collect the DFT energies as a dictionary:
        dft_energies = get_dft_energies(os.path.abspath(ilogfile))
        if not dft_energies:
            print("No dft energies found for index: ", logindex)
            failed_file_indices.append(logindex)
            continue
        # Get the smile string using the xyz2mol module.
        xyz_smiles = get_smiles_from_xyz(ixyzfile)
        if not xyz_smiles:
            print("Xyz not converted to file for index: ", logindex)
            update_failed_indices(failed_output, logindex)
            continue
        all_key_vals["index"] = logindex
        all_key_vals["smiles"] = xyz_smiles
        all_key_vals.update(dft_energies)
        save_to_output(output_csv, write_mode, all_key_vals)
        if i%10000 == 0:
            print("Number of DFT energies collected: ", (i+1))
    update_failed_indices(failed_output, failed_file_indices)
    return


if __name__ == "__main__":
    main()
