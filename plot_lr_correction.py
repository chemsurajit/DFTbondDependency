import argparse
import os
import glob
import pandas as pd


all_bonds = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
              'C_d_O', 'C_t_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N',
              'C_N_A', 'C_s_F', 'O_s_O', 'O_d_O', 'O_O_A', 'O_s_H',
              'O_s_N', 'O_d_N', 'O_t_N', 'O_N_A', 'O_s_F', 'N_s_N',
              'N_d_N', 'N_t_N', 'N_N_A', 'N_s_H', 'N_s_F', 'F_s_H']


def get_csv_files(csv_dir):
    """This function returns a list of csv files with matching pattern Reaction_*.csv."""
    files = []
    cwd = os.getcwd()
    os.chdir(csv_dir)
    for ifile in glob.glob('Reactions_*.csv'):
        files.append(os.path.abspath(ifile))
    os.chdir(cwd)
    return files


def get_energy_data(csv_files, lr_coeff_csv, dft_functional, bonds_ceffs = None):
    """Returns:
               dE, dE_corrected
    """
    dE = []
    dE_corrected = []
    chunksize = 100000
    print("bonds list: ", bonds_ceffs)
    print("Starting to read csv files")
    for csvs in csv_files:
        nchunk =0
        print("Reading csv file: ", csvs)
        for chunk in pd.read_csv(csvs, chunksize=chunksize):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            # remove all the reactions where there is no bond changes:
            df = df.loc[(df[all_bonds].abs().sum(axis=1) != 0)]
            dE_np = df.loc[:, ("G4MP2")] - df.loc[:, (dft_functional.upper())].to_numpy()
            print(dE_np)
            break
        break
    return 1, 2


def save_to_csv(dE, dE_corrected, output=None):
    return


def main():
    parser = argparse.ArgumentParser("Program to plot the correction from the result of linear regression.")
    parser.add_argument('-d', '--csv_directory', help="Directory where the Reactions_n.csv files are present",
                        required=True)
    parser.add_argument('-c', '--lr_coeff', help="CSV file that contains the coefficients from LR. Format: bonds,coeffs", required=True)
    parser.add_argument('-f', '--dft_functional', help="Name of the csv function for which the LR was performed.",
                        required=True)
    parser.add_argument('-n', '--name', help="Name of this run. Output files will be generated based on this name.",
                        required=True)
    args = parser.parse_args()
    output = args.name

    csv_files = get_csv_files(args.csv_directory)
    bonds_lr_pd = pd.read_csv(args.lr_coeff, index_col=False)
    dE, dE_corrected = get_energy_data(csv_files, args.lr_coeff, args.dft_functional, bonds_coeffs=dict(zip(bonds_lr_pd["bonds"], bonds_lr_pd["Coefficients"])))
    save_to_csv(dE, dE_corrected, output=output+".csv")
    pass


if __name__ == "__main__":
    main()
