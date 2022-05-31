import os, glob, sys
import argparse
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant


#bonds_list = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
#              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
#              'O_s_H', 'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
#              'N_N_A', 'N_s_H']

#remove the bonds with H
#bonds_list = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_O',
#              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
#              'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
#              'N_N_A']


bonds_list = ['C_d_C', 'C_t_C', 'C_C_A', 'C_s_O',
              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
              'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
              'N_N_A']

def detect_correlation(csv_files, output, frac=0.1):
    val_x = []
    val_y = []
    df_collect = []
    for csvf in csv_files:
        nchunk = 0
        print("Reading csv file: ", csvf)
        for chunk in pd.read_csv(csvf, chunksize=100000):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            df = df.loc[(df[bonds_list].abs().sum(axis=1) != 0)]
            if frac == 1.0:
                df_collect.append(df)
            else:
                df_collect.append(df.sample(frac=frac))
            if (nchunk%10) == 0:
                print("Nchunk:", nchunk)
    total_df = pd.concat(df_collect, ignore_index=True)
    print("Dataset reading complete. Proceeding for correlation check with ", frac, " of data")
    print("ndf row: ", total_df.shape[0])
    X = add_constant(total_df[bonds_list])
    # to return the result to a dataframe:
    vif_data = pd.DataFrame()
    vif_data["feature"] = X.columns
    vif_data["VIF"] = pd.Series(variance_inflation_factor(X.values, i) for i in range(X.shape[1]))
    print(vif_data)
    vif_data.to_csv(output)
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
    parser.add_argument('-f', '--frac_data', help="Fraction of data to be used to calculate the correlation. Value allowed: 0<frac<1.0.", default=0.5)

    args = parser.parse_args()
    csv_dir = args.data_dir
    csv_files = get_csv_files(csv_dir)
    frac = float(args.frac_data)
    output = "correlation.csv"
    if not output.endswith(".csv"):
        print("Only csv files allowed as output data file")
        parser.print_help()
        sys.exit()
    detect_correlation(csv_files, output, frac=frac)
    return


if __name__ == "__main__":
    main()
