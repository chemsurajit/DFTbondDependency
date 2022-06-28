import logging
import sys, os, glob
import pandas as pd
import numpy as np
import argparse
from sklearn.linear_model import LinearRegression
import statsmodels.api as sma
from scipy import stats

"""
Script to perform the linear regression.
"""

bonds_list = ['C_d_C', 'C_t_C', 'C_C_A', 'C_s_O',
              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
              'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
              'N_N_A']


def compute_statmodel_linear_regression(csv_files, dft_functional):
    """
    This function is for performing the linear regression using the statsmodel package. It is not used.
    :param csv_files: csv file containing reactions data
    :param dft_functional: Name of the DFT functional
    :return: None
    """
    val_x = []
    val_y = []
    for csvf in csv_files:
        nchunk = 0
        logging.info("Reading csv file: %s" % csvf)
        for chunk in pd.read_csv(csvf, chunksize=100000):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            # remove all the reactions where there is no bond changes:
            df = df.loc[(df[bonds_list].abs().sum(axis=1) != 0)]
            df["dE"] = df.loc[:, ("G4MP2")] - df.loc[:, (dft_functional.upper())]
            _val_x = df[bonds_list].values.tolist()
            if len(_val_x) > 0:
                val_x.extend(_val_x)
                val_y.extend(df["dE"].tolist())
            if nchunk%10 == 0:
                logging.info("Nchunk: %d" % nchunk)
        logging.info("Done reading file: %s" % csvf)
    if (len(val_x) > 0) and (len(val_y) > 0):
        val_x = np.array(val_x)
        val_y = np.array(val_y)
        val_x2 = sma.add_constant(val_x)
        est = sma.OLS(val_y, val_x2)
        est2 = est.fit()
        print(est2.summary())
    return


def compute_sklearn_linear_regression(csv_files, dft_functional, output):
    """
    Function to perform linear regression using sklearn.
    :param csv_files: csv file containing reaction data
    :param dft_functional: Name of the DFT functional
    :param output: csv filename to print out the result
    :return: None
    """
    val_x = []
    val_y = []
    for csvf in csv_files:
        nchunk = 0
        logging.info("Reading csv file: %s" % csvf)
        for chunk in pd.read_csv(csvf, chunksize=100000):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            # remove all the reactions where there is no bond changes:
            df = df.loc[(df[bonds_list].abs().sum(axis=1) != 0)]
            # This is the error between G4MP2 rean energies and the DFT rean energies.
            df["dE"] = df.loc[:, ("G4MP2")] - df.loc[:, (dft_functional.upper())]
            _val_x = df[bonds_list].values.tolist()
            if len(_val_x) > 0:
                val_x.extend(_val_x)
                val_y.extend(df["dE"].tolist())
            if nchunk%10 == 0:
                logging.info("Nchunk: %d" % nchunk)
        logging.info("Done reading file: %s" % csvf)
    if (len(val_x) > 0) and (len(val_y) > 0):
        val_x = np.array(val_x)
        val_y = np.array(val_y)
        logging.debug("Max and min of x: %f, %f" % (np.amax(val_x), np.amin(val_x)))
        logging.debug("Max and min of y: %f, %f" % (np.amax(val_y), np.amin(val_y)))
        reg = LinearRegression(fit_intercept=False).fit(val_x, val_y)
        logging.info("Regression score: %f " % reg.score(val_x, val_y))
        logging.info("The intercept: %f" % reg.intercept_)
        list_coef = reg.coef_.tolist()
        logging.info("bonds \t coefficient")
        logging.info("\n".join("{}\t{}".format(x, y) for x, y in zip(bonds_list, list_coef)))
        #
        # Calculation of P values.
        # source: https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression
        params = np.append(reg.intercept_, reg.coef_)
        predictions = reg.predict(val_x)
        newX = np.append(np.ones((len(val_x),1)), val_x, axis=1)
        MSE = (sum((val_y-predictions)**2))/(len(newX)-len(newX[0]))
        var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
        sd_b = np.sqrt(var_b)
        ts_b = params/sd_b
        p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b]
        lr_df = pd.DataFrame()
        logging.debug("Len bonds: %d, len coeffs: %d, len se: %d, len t: %d, len p: %d" %
                      (len(bonds_list), len(params), len(sd_b), len(ts_b), len(p_values)))
        lr_df["bonds"], lr_df["coefficients"], lr_df["standard_errors"], lr_df["t_values"], lr_df["probabilities"] = [bonds_list, params[1:], sd_b[1:], ts_b[1:], p_values[1:]]
        lr_df.to_csv(output, index=False)
    else:
        logging.info("Length of the x and y not greater than zero.")
        logging.info("Length bonds: %d" % len(val_x))
        logging.info("Length errors: %d" % len(val_y))
    logging.info("Linear regression done.")
    return


def get_csv_files(csv_dir, match=None):
    """
    This function returns a list of csv files with matching pattern Reaction_*.csv.
    """
    files = []
    cwd = os.getcwd()
    os.chdir(csv_dir)
    if match is None:
        match = "*.csv"
    for ifile in glob.glob(match):
        files.append(os.path.abspath(ifile))
    os.chdir(cwd)
    return files


def get_arguments():
    """
    Function to process arguements.
    """
    parser = argparse.ArgumentParser(
        "Script to perform linear regression for different DFT functionals."
    )
    parser.add_argument(
        "-data_dir", "--data_dir",
        type=str,
        help="Path where the Reactions_n.csv files are kept. \
                        default, current directory.",
        required=False,
        default="./"
    )

    parser.add_argument(
        "-dft_functional", "--dft_functional",
        type=str,
        help="Name of the dft functional \
                        from which to calculate errors in reaction energies. Required only if --limit is used.",
        required=True,
    )

    parser.add_argument(
        "-output", "--output",
        type=str,
        help="Output csv file for the coefficients and other values.",
        required=False,
        default="output.csv"
    )
    parser.add_argument(
        "-logging", "--logging",
        type=str,
        required=False,
        default="info",
        choices=["debug", "info", "warning", "error", "critical"],
        help="Provide logging level. Default is warning."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_arguments()
    log_level = args.logging.upper()
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=log_level,
        datefmt="%H:%M:%S",
    )
    csv_dir = args.data_dir
    match = "Reactions_*.csv"
    logging.info("Getting all the data csv files from %s" % csv_dir)
    csv_files = get_csv_files(csv_dir, match=match)
    dft_functional = args.dft_functional
    logging.debug("DFT method for the LR: %s" % dft_functional)
    output = args.output
    if not output.endswith(".csv"):
        logging.error("Only csv files allowed as output data file")
        sys.exit()
    logging.info("Now performing the linear regression.")
    compute_sklearn_linear_regression(csv_files, dft_functional, output)
    #compute_statmodel_linear_regression(csv_files, dft_functional, output)
