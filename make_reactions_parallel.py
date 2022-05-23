import os
import argparse
import sys
import time
import pandas as pd
import csv
import numpy as np
import concurrent.futures as confut
import logging

bonds_list = [
            'C_s_C', 'C_d_C', 'C_t_C', 'C_C_A', 'C_s_H', 'C_s_O',
            'C_d_O', 'C_t_O', 'C_O_A', 'C_s_N', 'C_d_N', 'C_t_N',
            'C_N_A', 'C_s_F', 'O_s_O', 'O_d_O', 'O_O_A', 'O_s_H',
            'O_s_N', 'O_d_N', 'O_t_N', 'O_N_A', 'O_s_F', 'N_s_N',
            'N_d_N', 'N_t_N', 'N_N_A', 'N_s_H', 'N_s_F', 'F_s_H'
]



def get_required_mol_data(mol_data_file):
    """
    This function returns the molecules data as pandas dataframe.
    Only the index, smiles, chemformula, bonds, B3LYP-D, PBE, M06-2X columns will be returned.
    :param mol_data_file:
    :return: mold_data as pandas dataframe
    """
    other_columns = [
        "index", "smiles", "chemformula", "B3LYP-D", "PBE", "M06-2X"
    ]
    columns_to_read = bonds_list + other_columns
    mol_data_required = pd.read_csv(mol_data_file, usecols=columns_to_read,
                                    keep_default_na=False, na_values=np.nan)
    return mol_data_required.dropna()


def get_required_g4mp2_data(g4mp2_en_file):
    """
    This function will read a csv file containing G4MP2 energy values and index
    as two must present columns and return those two columns as dictionary.
    :param g4mp2_en_file:
    :return: G4MP2 dataframe as {index:energy}
    """
    pdata = pd.read_csv(g4mp2_en_file, usecols=["index","G4MP2"],
                        keep_default_na=False, na_values=np.nan,
                        ).dropna()
    return pdata


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Make csv files containing all the reactions data for the three DFT functional."
    )
    parser.add_argument(
        "-rid_csv", "--rid_csv",
        type=str,
        required=True,
        help="csv file containing reactant and product indices"
    )
    parser.add_argument(
        "-mol_data", "--mol_data",
        type=str,
        required=True,
        help="csv file containing all QM9 molecules data"
    )
    parser.add_argument(
        "-g4mp2_en", "--g4mp2_en",
        type=str,
        required=True,
        help="csv file containing index of molecules and G4MP2 energies."
    )
    parser.add_argument(
        "-nprocs", "--nprocs",
        type=int,
        required=False,
        default=1,
        help="Number of processors to be used."
    )
    parser.add_argument(
        "-out_dir", "--out_dir",
        type=str,
        required=False,
        default="../outputs",
        help="Directory where the Reactions_n.csv files will be saved."
    )
    parser.add_argument(
        '-log', '--log',
        type=str,
        required=False,
        default="warning",
        choices=["debug", "info", "warning", "error", "critical"],
        help="Provide logging level. Default is warning."
    )
    return parser.parse_args()


def process_reaction_data(rids_pd, outid, molecule_data_pd, g4mp2_en, outdir):
    """
    The main function where all the reactions will be computed.
    :param rids_pd, outid, molecules_pd, g4mp2_pd, outdir
    :return: Integer 0 upon completion.
    """
    logging.debug("Inside process_reaction_data function")
    bonds_ens_cols = bonds_list + ["PBE","B3LYP-D","M06-2X"]
    output_csv_file = os.path.join(outdir, "Reactions_" + str(outid) + ".csv")
    pid = os.getpid()
    ppid = os.getppid()
    logging.info("start index: %d, pid: %d" % (list(rids_pd.index.values)[0], pid))
    logging.info("  End index: %d, pid: %d" % (list(rids_pd.index.values)[-1], pid))
    logging.info("pid: %d, nreaction to be processed: %d" % (pid, rids_pd.shape[0]))
    start = time.time()
    logging.info("pid, ppid info: %s %s" % (pid, ppid))
    #do stuff here.
    counter = 0
    for rowid, row in rids_pd.iterrows():
        logging.debug("loopstart pid, rowid, row: %d, %d, %s" % (pid, rowid, row.to_string()))
        reactant_index = row.reactindex
        logging.debug("pid, reactant index: %d %d" % (pid, reactant_index))
        pdt_index = row.pdtindex
        logging.debug("pid, pdt index: %d %d" % (pid, pdt_index))
        try:
            reactant_row = molecule_data_pd.loc[molecule_data_pd['index'] == reactant_index]
        except Exception as ex:
            logging.warning("The following exception occurs: %s " % str(ex))
            logging.warning("No row in molecule_data_pd for pid, reactant index: %d %d" % (pid, reactant_index))
            continue
        logging.debug("pid, react row: %d %s" % (pid, reactant_row.to_string()))
        try:
            pdt_row = molecule_data_pd.loc[molecule_data_pd['index'] == pdt_index]
        except Exception as ex:
            logging.warning("The following exception occurs: %s " % str(ex))
            logging.warning("No row in molecule_data_pd for pid, pdt index: %d %d" % (pid, pdt_index))
            continue
        logging.debug("pid, pdt row: %d %s" % (pid, pdt_row.to_string()))
        react_smi = reactant_row["smiles"].values[0]
        logging.debug("pid, react_smi %d %s" % (pid, react_smi))
        pdt_smi = pdt_row["smiles"].values[0]
        logging.debug("pid, pdt_smi: %d %s" % (pid, pdt_smi))
        reaction_prop_diff = pdt_row[bonds_ens_cols] - reactant_row[bonds_ens_cols].values
        logging.debug("pid, reaction_prop_diff1: %d %s" % (pid, reaction_prop_diff))
        reaction_prop_diff["react_smi"], reaction_prop_diff["pdt_smi"] = [react_smi, pdt_smi]
        logging.debug("pid, reaction_prop_diff2: %d, %s" % (pid, reaction_prop_diff.to_string()))
        try:
            react_g4mp2_en = g4mp2_en.loc[g4mp2_en["index"] == reactant_index, "G4MP2"].values[0]
        except Exception as ex:
            logging.warning("The following exception occurs: %s " % str(ex))
            logging.warning("No G4MP2 energy found for pid, reactant index: %d %d" % (pid, reactant_index))
            continue
        try:
            pdt_g4mp2_en = g4mp2_en.loc[g4mp2_en["index"] == pdt_index, "G4MP2"].values[0]
        except Exception as ex:
            logging.warning("The following exception occurs: %s " % str(ex))
            logging.warning("No G4MP2 energy found for pid, pdtindex: %d %d" % (pid, pdt_index))
            continue
        reaction_prop_diff["G4MP2"] = pdt_g4mp2_en - react_g4mp2_en
        logging.debug("pid, reaction_prop_diff3: %d, %s" %(pid, reaction_prop_diff.to_string()))
        reaction_prop_diff["chemformula"] = reactant_row["chemformula"].values[0]
        logging.debug("pid, reaction_prop_diff4: %d, %s" % (pid, reaction_prop_diff.to_string()))
        reaction_prop_diff["reactindex"], reaction_prop_diff["pdtindex"] = [reactant_index, pdt_index]
        logging.debug("pid, reaction_prop_diff5: %d, %s" % (pid, reaction_prop_diff.to_string()))
        logging.debug("loopend pid, rowid: %d, %d" % (pid, rowid))
        if counter == 0:
            logging.info("New csv file will be created file name: %s, pid: %d" % (output_csv_file, pid))
            logging.debug("csv, pid, rowid in if: %s, %d, %d" %(output_csv_file, pid, rowid))
            logging.debug(reaction_prop_diff.keys())
            reaction_prop_diff.to_csv(output_csv_file,
                                      mode='w', index=False,
                                      quoting=csv.QUOTE_MINIMAL,
                                      sep=",")
        else:
            logging.debug("Else, append to csv pid, rowid: %s, %d, %d" % (output_csv_file, pid, rowid))
            reaction_prop_diff.to_csv(output_csv_file,
                                      mode='a', index=False,
                                      quoting=csv.QUOTE_MINIMAL,
                                      header=False, sep=",")
            logging.debug("file %s updated" % output_csv_file)
        if (counter+1) % 100000 == 0:
            logging.info("Converted: %d reactions to %s with pid %d" % (counter, output_csv_file, pid))
        counter += 1
    stop = time.time()
    completed_in = round((stop-start)/3600.0, 2)
    logging.info("Loop with pid %d completed in: %s hr" % (pid, completed_in))
    return


if __name__ == "__main__":
    # multi threading credit: https://betterprogramming.pub/pandas-how-to-process-a-dataframe-in-parallel-make-pandas-lightning-fast-669978cf5356
    args = get_arguments()
    log_level = args.log.upper()
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=log_level,
        datefmt="%H:%M:%S",
    )
    outdir = os.path.abspath(args.out_dir)
    nprocs = args.nprocs
    logging.info("Number of processors chose: %s" % nprocs)
    molecule_data_pd = get_required_mol_data(args.mol_data)
    g4mp2_en = get_required_g4mp2_data(args.g4mp2_en)
    total_reaction_ids_pd = pd.read_csv(args.rid_csv,
                                        usecols=["reactindex","pdtindex"],
                                        keep_default_na=False, na_values=np.nan
                                        ).dropna()
    splitted_rid_pd = np.array_split(total_reaction_ids_pd, nprocs)
    # setup related to multiprocessing
    main_func_results = []
    start = time.time()
    logging.info("Starting parallel run.")
    with confut.ProcessPoolExecutor(max_workers=nprocs) as executor:
        results = [executor.submit(process_reaction_data, rid_pd, npd, molecule_data_pd, g4mp2_en, outdir) for npd, rid_pd in enumerate(splitted_rid_pd)]
        for result in confut.as_completed(results):
            try:
                main_func_results.append(result.result())
            except Exception as ex:
                # taken from: https://www.codegrepper.com/code-examples/python/python+exception+with+line+number
                logging.error("Exception occurs: %s " % str(ex))
                ex_type, ex_obj, ex_trace = sys.exc_info()
                line_number = ex_trace.tb_lineno
                logging.error("ERROR: %s. LINE NO: %s" % (str(ex), line_number))
                pass
    end = time.time()
    logging.info("JOB COMPLETED.")
    logging.info("PPID %s Completed in %s hr" % (os.getpid(), round((end-start)/3600,2)))
