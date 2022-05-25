import json
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
    parser.add_argument(
        '-i', '--indices',
        type=str,
        required=True,
        help="This file should contains the indices of the datafile to be read. Format should be in json."
    )
    return parser.parse_args()


def process_reaction_data_column(rids_pd, coreno, nodeno, molecule_data_pd, g4mp2_en, outdir):
    logging.debug("Inside process_reaction_data function")
    bonds_ens_cols = bonds_list + ["PBE", "B3LYP-D", "M06-2X"]
    output_csv_file = os.path.join(outdir, "Reactions_" + str(nodeno) + "_core_" + str(coreno) + ".csv")
    pid = os.getpid()
    ppid = os.getppid()
    logging.info("start index: %d, pid: %d" % (list(rids_pd.index.values)[0], pid))
    logging.info("  End index: %d, pid: %d" % (list(rids_pd.index.values)[-1], pid))
    logging.info("pid: %d, nreaction to be processed: %d" % (pid, rids_pd.shape[0]))
    start = time.time()
    logging.info("pid, ppid info: %s %s" % (pid, ppid))
    nchunks = 0
    chunksize = rids_pd.shape[0] // 10000
    logging.debug("Start processing data for pid: %d in for loop" % pid)
    for chunk in np.array_split(rids_pd, chunksize):
        react_indices = chunk["reactindex"].to_list()
        pdt_indices = chunk["pdtindex"].to_list()
        logging.info("length of reactant and product indices %d %d. pid: %d" % (len(react_indices), len(pdt_indices), pid))
        react_g4mp2_ens = g4mp2_en.loc[react_indices, "G4MP2"]
        pdt_g4mp2_ens = g4mp2_en.loc[pdt_indices, "G4MP2"]
        print(react_g4mp2_ens.head())
        logging.info("No of react g4mp2 energy rows: %d. pid: %d" % (react_g4mp2_ens.shape[0], pid))
        logging.info("No of pdt g4mp2 energy rows: %d. pid: %d" % (pdt_g4mp2_ens.shape[0], pid))
        mol_react_data = molecule_data_pd.loc[react_indices, bonds_ens_cols]
        mol_pdt_data = molecule_data_pd.loc[pdt_indices, bonds_ens_cols]
        logging.debug("number of rows in reactant df: %d. pid: %d" % (mol_react_data.shape[0], pid))
        logging.debug("number of rows in pdt df: %d. pid: %d" % (mol_pdt_data.shape[0], pid))
        logging.info("Now dataframe wise subtraction will be performed pid: %d" % pid)
        reaction_result_pd = mol_pdt_data[bonds_ens_cols].subtract(mol_react_data[bonds_ens_cols])
        logging.debug("done subtracting dataframe. pid: %d" % pid)
        reaction_result_pd["G4MP2"] = pdt_g4mp2_ens["G4MP2"] - react_g4mp2_ens["G4MP2"]
        logging.debug("Done updating G4MP2 energy. pid: %d" % pid)
        reaction_result_pd["chemformula"] = mol_react_data["chemformula"]
        logging.debug("Updated chemformula to reaction data. pid: %d" % pid)
        reaction_result_pd["react_smi"] = mol_react_data["smiles"]
        logging.debug("Updated react_smi to reaction data. pid: %d" % pid)
        reaction_result_pd["pdt_smi"] = mol_pdt_data["smiles"]
        logging.debug("Updated pdt_smi to reaction data. pid: %d" % pid)
        reaction_result_pd["reactindex"] = react_indices
        logging.debug("Updated reactindex to reaction data. pid: %d" % pid)
        reaction_result_pd["pdtindex"] = pdt_indices
        logging.debug("Updated pdtindex to reaction data. pid: %d" % pid)
        if nchunks == 0:
            logging.info("New csv file will be created file name: %s, pid: %d" % (output_csv_file, pid))
            reaction_result_pd.to_csv(output_csv_file,
                                      mode="w", index=False,
                                      quoting=csv.QUOTE_MINIMAL,
                                      sep=",")
        else:
            logging.info("Will update the datachunk to csv file: %s" % output_csv_file)
            reaction_result_pd.to_csv(output_csv_file,
                                      mode='a', index=False,
                                      quoting=csv.QUOTE_MINIMAL,
                                      header=False, sep=",")
        nchunks += 1
        if nchunks%10 == 0:
            logging.info("Converted: %d reactions to %s with pid %d" % (nchunks, output_csv_file, pid))
    stop = time.time()
    completed_in = round((stop-start)/3600.0, 2)
    logging.info("Loop with pid %d completed in: %s hr" % (pid, completed_in))
    return


def process_reaction_data(rids_pd, coreno, nodeno, molecule_data_pd, g4mp2_en, outdir):
    """
    The main function where all the reactions will be computed.
    :param rids_pd, outid, molecules_pd, g4mp2_pd, outdir
    :return: Integer 0 upon completion.
    """
    logging.debug("Inside process_reaction_data function")
    bonds_ens_cols = bonds_list + ["PBE", "B3LYP-D", "M06-2X"]
    output_csv_file = os.path.join(outdir, "Reactions_" + str(nodeno) + "_core_" + str(coreno) + ".csv")
    pid = os.getpid()
    ppid = os.getppid()
    logging.info("start index: %d, pid: %d" % (list(rids_pd.index.values)[0], pid))
    logging.info("  End index: %d, pid: %d" % (list(rids_pd.index.values)[-1], pid))
    logging.info("pid: %d, nreaction to be processed: %d" % (pid, rids_pd.shape[0]))
    start = time.time()
    logging.info("pid, ppid info: %s %s" % (pid, ppid))
    #do stuff here.
    counter = 0
    chunk_tocsv = []
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
            counter += 1
            continue
        # the dataframe will only be appended from here to skip the first row.
        chunk_tocsv.append(reaction_prop_diff)
        #else:
        #    logging.debug("Else, append to csv pid, rowid: %s, %d, %d" % (output_csv_file, pid, rowid))
        #    reaction_prop_diff.to_csv(output_csv_file,
        #                              mode='a', index=False,
        #                              quoting=csv.QUOTE_MINIMAL,
        #                              header=False, sep=",")
        #    logging.debug("file %s updated" % output_csv_file)
        if (counter+1) % 50000 == 0:
            logging.info("Converted: %d reactions to %s with pid %d" % (counter, output_csv_file, pid))
            logging.info("Will update the datachunk to csv file: %s" % output_csv_file)
            pd.concat(chunk_tocsv).to_csv(output_csv_file,
                                      mode='a', index=False,
                                      quoting=csv.QUOTE_MINIMAL,
                                      header=False, sep=",")
            chunk_tocsv = []
        counter += 1
    # Now write the remaining data at the end to the csv file
    logging.info("Will update the reminder of the datachunk outside for loop to %s. pid: %d" % (output_csv_file, pid))
    pd.concat(chunk_tocsv).to_csv(output_csv_file,
                                  mode='a', index=False,
                                  quoting=csv.QUOTE_MINIMAL,
                                  header=False, sep=",")
    stop = time.time()
    completed_in = round((stop-start)/3600.0, 2)
    logging.info("Loop with pid %d completed in: %s hr" % (pid, completed_in))
    return


def get_big_chunk_rid_pd(reaction_csv_file, indices_file):
    with open(indices_file) as fp:
        json_dict = json.load(fp)
    index_list = json_dict["indices"]
    del json_dict
    df = pd.read_csv(reaction_csv_file,
                     keep_default_na=False, na_values=np.nan
                     ).dropna()
    selected_df = df.iloc[index_list]
    del df
    return selected_df


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
    node_no = args.indices.split(".")[0]
    logging.info("Number of processors chose: %s" % nprocs)
    # load the indices of the reaction as total data frame
    logging.info("loading %s ..." % args.rid_csv)
    bigchunk_rid_pd = get_big_chunk_rid_pd(args.rid_csv, args.indices)
    # load the molecular bond and DFT energy data
    logging.info("loading %s ..." % args.mol_data)
    molecule_data_pd = get_required_mol_data(args.mol_data)
    # load the g4mp2 molecular data
    logging.info("loading %s ..." % args.g4mp2_en)
    g4mp2_en = get_required_g4mp2_data(args.g4mp2_en)
    # modify the molecule_data_pd according to the indices of the G4MP2 energy.
    logging.info("The molecule data will be taken for the indices of the G4MP2 database.")
    g4mp2_indices = g4mp2_en["index"].to_list()
    modified_mol_data_pd = molecule_data_pd.loc[molecule_data_pd["index"].isin(g4mp2_indices)]
    del molecule_data_pd
    #total_reaction_ids_pd = pd.read_csv(args.rid_csv,
    #                                    usecols=["reactindex","pdtindex"],
    #                                    keep_default_na=False, na_values=np.nan
    #                                    ).dropna()
    splitted_rid_pd = np.array_split(bigchunk_rid_pd, nprocs)
    # free the memory where the larger reaction.csv data is loaded.
    del bigchunk_rid_pd
    # setup related to multiprocessing
    main_func_results = []
    start = time.time()
    logging.info("Starting parallel run.")
    with confut.ProcessPoolExecutor(max_workers=nprocs) as executor:
        results = [executor.submit(process_reaction_data, rid_pd, coreno, node_no, modified_mol_data_pd, g4mp2_en, outdir) for coreno, rid_pd in enumerate(splitted_rid_pd)]
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
