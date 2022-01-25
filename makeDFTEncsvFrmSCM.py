import os, sys
import fnmatch
import csv
import pandas as pd
import argparse


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
            logfiles.append(files)
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
    for file in os.listdir(directory):
        # Take the xyz files which has 'dsgdb9nsd_' in their name
        if fnmatch.fnmatch(file, 'dsgdb9nsd_*.xyz'):
            xyz_files.append(file)
    # Check if the list xyz_files is empty. If it is empty, then program exits.
    if not xyz_files:
        print("No files with extension .xyz found in directory: ", directory)
        print("Program exit now.")
        sys.exit()
    return sorted(xyz_files)


def main():
    parser = argparse.ArgumentParser("Program to create csv files from the logfiles of SCM output.")
    parser.add_argument('-d', '--dir_logfiles',
                        help="Location of the directory from where the logfiles will be read.",
                        required=False, default='./')
    parser.add_argument('-xd', '--xyzd', help="Directory where the xyz files are present.",
                        required=True)
    parser.add_argument('-o', '--output_csv', help="Name of the output csv file for the energies.",
                        required=False, default='qm9_dft_en.csv')
    args = parser.parse_args()
    xyz_dir = args.xyzd
    output_csv = args.output_csv
    log_dir = args.dir_logfiles
    print("XYZ dir: ", xyz_dir)
    print("output csv: ", output_csv)
    print("log directory: ", log_dir)
    logfiles = get_all_logfiles(log_dir)
    print("Number of logfiles: ", len(logfiles))
    xyzfiles = get_xyz_files(xyz_dir)
    print("Number of xyzfiles: ", len(xyzfiles))
    return


if __name__ == "__main__":
    main()
