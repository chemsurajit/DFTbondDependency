import csv
import pandas as pd
import argparse


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
    return


if __name__ == "__main__":
    main()
