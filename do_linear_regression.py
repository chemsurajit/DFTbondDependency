import sys, os, glob
import pandas as pd
import numpy as np
import argparse
from sklearn.linear_model import LinearRegression
from scipy import stats


bonds_list = ['C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
              'C_d_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N', 'C_N_A',
              'O_s_H', 'O_s_N', 'O_d_N', 'O_N_A', 'N_s_N', 'N_d_N',
              'N_N_A', 'N_s_H']


def compute_sklearn_linear_regression(csv_files, dft_functional, output):
    val_x = []
    val_y = []
    for csvf in csv_files:
        nchunk = 0
        print("Reading csv file: ", csvf)
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
                print("Nchunk: ", nchunk)
        print("Done reading file: ", csvf)
    if (len(val_x) > 0) and (len(val_y) > 0):
        val_x = np.array(val_x)
        val_y = np.array(val_y)
        print("Max and min of x: ", np.amax(val_x), np.amin(val_x))
        print("Max and min of y: ", np.amax(val_y), np.amin(val_y))
        reg = LinearRegression().fit(val_x, val_y)
        print("Regression score: ", reg.score(val_x, val_y))
        print("The intercept: ", reg.intercept_)
        #
        # Calculation of P values.
        # source: https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression
        params = np.append(reg.intercept_, reg.coef_)
        predictions = reg.predict(val_x)
        newX = np.append(np.ones((len(val_x),1)), val_x, axis=1)
        MSE = (sum((val_y-predictions)**2))/(len(newX)-len(newX[0]))
        var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
        sd_b = np.sqrt(var_b)
        ts_b = params/ sd_b
        p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b]
        myDF3 = pd.DataFrame()
        myDF3["Coefficients"],myDF3["Standard Errors"],myDF3["t values"],myDF3["Probabilities"] = [params,sd_b,ts_b,p_values]
        myDF3.to_csv(output)
    else:
        print("Length of the x and y not greater than zero.")
        print("Length x: ", len(val_x))
        print("Length y: ", len(val_y))
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
