import argparse
import glob
import sys
import os
import pandas as pd
from collections import Counter


bonds_list = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
              'C_d_O', 'C_t_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N',
              'C_N_A', 'C_s_F', 'O_s_O', 'O_d_O', 'O_O_A', 'O_s_H',
              'O_s_N', 'O_d_N', 'O_t_N', 'O_N_A', 'O_s_F', 'N_s_N',
              'N_d_N', 'N_t_N', 'N_N_A', 'N_s_H', 'N_s_F', 'F_s_H']


dft_functionals = ['G4MP2', 'KCIS-MODIFIED', 'KCIS-ORIGINAL', 'PKZB',
                   'VS98', 'LDA(VWN)', 'PW91', 'BLYP', 'BP', 'PBE', 'RPBE',
                   'REVPBE', 'OLYP', 'FT97', 'BLAP3', 'HCTH/93', 'HCTH/120',
                   'HCTH/147', 'HCTH/407', 'BMTAU1', 'BOP', 'PKZBX-KCISCOR',
                   'VS98-X(XC)', 'VS98-X-ONLY', 'BECKE00', 'BECKE00X(XC)',
                   'BECKE00-X-ONLY', 'BECKE88X+BR89C', 'OLAP3', 'TPSS', 'MPBE',
                   'OPBE', 'OPERDEW', 'MPBEKCIS', 'MPW', 'TAU-HCTH', 'XLYP', 'KT1',
                   'KT2', 'M06-L', 'BLYP-D', 'BP86-D', 'PBE-D', 'TPSS-D', 'B97-D',
                   'REVTPSS', 'PBESOL', 'RGE2', 'SSB-D', 'MVS', 'MVSX', 'T-MGGA',
                   'TPSSH', 'B3LYP(VWN5)', 'O3LYP(VWN5)', 'KMLYP(VWN5)', 'PBE0',
                   'B3LYP*(VWN5)', 'BHANDH', 'BHANDHLYP', 'B97', 'B97-1', 'B97-2',
                   'MPBE0KCIS', 'MPBE1KCIS', 'B1LYP(VWN5)', 'B1PW91(VWN5)', 'MPW1PW',
                   'MPW1K', 'TAU-HCTH-HYBRID', 'X3LYP(VWN5)', 'OPBE0', 'M05', 'M05-2X',
                   'M06', 'M06-2X', 'B3LYP-D']


def count_bond_change(csv_files, output=None):
    print(os.getcwd())
    total_sumed_counter = Counter({})
    for csvf in csv_files:
        nchunk = 0
        print("Reading csv file: ", csvf)
        for chunk in pd.read_csv(csvf, chunksize=100000):
            df = chunk.dropna(axis=0, how='any')
            abs_sum_counter = Counter(df[bonds_list].abs().sum(axis=0).to_dict())
            total_sumed_counter.update(abs_sum_counter)
            nchunk += 1
            if nchunk % 5 == 0:
                print("Nchunk: ", nchunk)
        print("Done reading csvfile: ", csvf)
    # save to csv file with format:
    # bonds,count
    with open(output, 'w') as fp:
        fp.write("bonds,count\n")
        for key, value in total_sumed_counter.items():
            fp.write('{},{}\n'.format(key, int(value)))
    return


def plot_zero_bond(csv_files, output=None, dft_func=None):
    print("CWD: ", os.getcwd())
    mode = 'w'
    if dft_func is None:
        print("No dft functional is supplied to calculate DFT errors.")
    else:
        selected_columns = ["reactionindex", "G4MP2", dft_func.upper(), "dE"]
        for csvf in csv_files:
            nchunk = 0
            print("Reading csv file: ", csvf)
            for chunk in pd.read_csv(csvf, chunksize=100000):
                nchunk += 1
                df = chunk.dropna(axis=0, how='any')
                #selected_columns_df = df[["reactionindex", dft_func.upper(), "G4MP2"]]
                #selected_columns_df["dE"] = df.loc[:, ("G4MP2")] - df.loc[:, (dft_func.upper())]
                df["dE"] = df.loc[:, ("G4MP2")] - df.loc[:, (dft_func.upper())]
                if mode == 'w':
                    df.to_csv(output, mode=mode, index=False, columns = selected_columns)
                    #selected_columns_df.to_csv(output, mode=mode, index=False)
                elif mode == 'a':
                    df.to_csv(output, mode=mode, index=False, header=False, columns=selected_columns)
                    #selected_columns_df.to_csv(output, mode=mode, index=False, header=False)
                mode = 'a'
                if nchunk % 10 == 0:
                    print("Nchunk: ", nchunk)
            print("Reading complete for file: ", csvf)
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
    parser = argparse.ArgumentParser("Analysis script for the reaction data.")
    parser.add_argument("-b", "--bond_change", help="To plot the number of bond change for \
                        all the reactions for each type of bonds", action='store_true')

    parser.add_argument("-d", "--data_dir", help="Path where the Reactions_n.csv files are kept. \
                        default, current directory.", default="./")

    parser.add_argument("-z", "--zero_de", help="To create a csv file for dE(DFT error) \
                        for reactions with zero bond change.", action='store_true')

    parser.add_argument('-f', '--dft_functional', required='--zero_de' in sys.argv, help="Name of the dft functional \
                        from which to calculate errors in reaction energies. Required only if --limit is used.")

    parser.add_argument('-o', '--output', required=True, help="The output to be printed. It needs to be image format file.")

    args = parser.parse_args()
    csv_dir = args.data_dir
    csv_files = get_csv_files(csv_dir)
    #
    # calling of functions
    if args.bond_change:
        print("Calculating bond change for all the reactions")
        output = args.output
        if not output.endswith(".csv"):
            print("Only csv file allowed as output data file.")
            sys.exit()
        count_bond_change(csv_files, output=output)
    if args.zero_de:
        print("Calculating DFT error for reactions with zero bond change.")
        output = args.output
        dft_functional = args.dft_functional
        if not output.endswith(".csv"):
            print("Only csv files allowed as output data file")
        plot_zero_bond(csv_files, output=output, dft_func=dft_functional)
    return


if __name__ == "__main__":
    main()
