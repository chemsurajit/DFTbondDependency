import logging
import os, glob, sys
import argparse
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant


"""
This script is to detect correlation between variables using variance inflation factor.
"""

# This list of bonds are hard coded as only these bonds are found to
# be changing in the QM9 dataset.
bonds_list = ['C_d_C', 'C_t_C', 'C_C_A', 'C_s_O',
              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N',
              'C_t_N', 'C_N_A', 'O_s_N', 'O_d_N',
              'O_N_A', 'N_s_N', 'N_d_N', 'N_N_A']

def detect_correlation(csv_files, output, frac=0.1):
    """
    The main function.
    :param csv_files: list of csv files containing reaction data
    :param output: csv filename to save the result
    :param frac: fraction of data to be used for the calculation
    :return: None
    """
    df_collect = []
    for csvf in csv_files:
        nchunk = 0
        logging.info("Reading csv file: %s" % csvf)
        for chunk in pd.read_csv(csvf, chunksize=100000):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            df = df[bonds_list]
            df = df.loc[(df[bonds_list].abs().sum(axis=1) != 0)]
            if frac == 1.0:
                df_collect.append(df)
            else:
                df_collect.append(df.sample(frac=frac))
            if (nchunk%10) == 0:
                logging.info("Nchunk:", nchunk)
    total_df = pd.concat(df_collect, ignore_index=True)
    logging.info("Dataset reading complete. Proceeding for correlation check with %f of the data." % frac)
    logging.info("Total number of reactions: %d" % total_df.shape[0])
    X = add_constant(total_df[bonds_list])
    # to return the result to a dataframe:
    vif_data = pd.DataFrame()
    vif_data["feature"] = X.columns
    logging.info("Performing VIF...")
    vif_data["VIF"] = pd.Series(variance_inflation_factor(X.values, i) for i in range(X.shape[1]))
    logging.info("VIF done. saving values to %s." % output)
    print(vif_data)
    vif_data.to_csv(output, index=False)
    return


def get_csv_files(csv_dir, match=None):
    """
    This function returns a list of csv files with matching pattern Reaction_*.csv.
    """
    files = []
    cwd = os.getcwd()
    # if no match provided, all csv files will be read.
    if match is None:
        match = "*.csv"
    os.chdir(csv_dir)
    for ifile in glob.glob(match):
        files.append(os.path.abspath(ifile))
    os.chdir(cwd)
    return files


def get_arguments():
    """
    function to process arguments.
    """
    parser = argparse.ArgumentParser("Script to perform correlation detection using VIF.")
    parser.add_argument(
        "-data_dir", "--data_dir",
        type=str,
        help="Path where the Reactions_n.csv files are kept. \
            default, current directory.",
        required=True,
        default="./"
    )
    parser.add_argument(
        "-frac_data", "--frac_data",
        type=float,
        help="Fraction of data to be used to calculate the correlation. Value allowed: 0<frac<1.0.",
        required=False,
        default=0.5
    )
    parser.add_argument(
        "-output", "--output",
        type=str,
        help="Output file to save the correlation values. CSV format.",
        required=False,
        default="correlation.csv"
    )
    return parser.parse_args()

def main():
    args = get_arguments()
    csv_dir = args.data_dir
    match = "Reactions_*.csv"
    csv_files = get_csv_files(csv_dir, match=match)
    frac = float(args.frac_data)
    output = "correlation.csv"
    if not output.endswith(".csv"):
        logging.error("Only csv files allowed as output data file")
        parser.print_help()
        sys.exit()
    detect_correlation(csv_files, output, frac=frac)
    return


if __name__ == "__main__":
    main()
