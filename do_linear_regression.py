import sys, os, glob
import pandas as pd
import numpy as np
import argparse
from sklearn.linear_model import LinearRegression


bonds_list = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
              'O_s_H', 'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
              'N_N_A', 'N_s_H']


def compute_sklearn_linear_regression(csv_files, dft_functional, output):
    for csvf in csv_files:
        nchunk = 0
        for chunk in pd.read_csv(csvf, chunksize=100000):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            # remove all the reactions where there is no bond changes:

            df["dE"] = df.loc[:, ("G4MP2")] - df.loc[:, (dft_functional.upper())]
    return


def get_csv_files(csv_dir):
    """This function returns a list of csv files with matching pattern Reaction_*.csv."""
    files = []
    cwd = os.getcwd()
    os.chdir(csv_dir)
    for ifile in glob.glob('Reactions_*.csv'):
        files.append(os.path.abspath(ifile))
    os.chdir(cwd)
    return files


def main():
    parser = argparse.ArgumentParser("Script to perform linear regression for different DFT functionals.")
    parser.add_argument('-d', '--data_dir', help="Path where the Reactions_n.csv files are kept. \
                        default, current directory.", default="./")

    parser.add_argument('-f', '--dft_functional', required=True, help="Name of the dft functional \
                        from which to calculate errors in reaction energies. Required only if --limit is used.")

    parser.add_argument('-o', '--output', help="Output csv file for the coefficients and other values.", default="output.csv")

    args = parser.parse_args()
    csv_dir = args.data_dir
    csv_files = get_csv_files(csv_dir)
    dft_functional = args.dft_functional

    output = args.output
    if not output.endswith(".csv"):
        print("Only csv files allowed as output data file")
        parser.print_help()
        sys.exit()

    compute_sklearn_linear_regression(csv_files, dft_functional, output)

    return


if __name__ == "__main__":
    main()
