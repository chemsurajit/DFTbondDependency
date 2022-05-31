from ase.formula import Formula as aseformula
import argparse
import pandas as pd
import csv
import logging
import time


def find_reactions(molecule_csv_file, output=""):
    """function to create csv files with reactant and product ids."""
    # reading from csv file containing molecules data but only 4 cols.
    mol_csv_data = pd.read_csv(molecule_csv_file, usecols=["index", "smiles", "chemformula"])
    indices = mol_csv_data["index"].tolist()
    formulas = mol_csv_data["chemformula"].tolist()
    nmols = len(indices)
    logging.info("Reaction search will start now.")
    with open(output, 'w') as fcsv:
        writer = csv.writer(fcsv, delimiter=",", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerow(['reactindex', 'pdtindex'])
        reaction_ids = []
        for i in range(nmols):
            reactid = indices[i]
            reactformula = formulas[i]
            # the following return dictionary of number of each element in a molecule
            react_elem_count = aseformula(reactformula).count()
            for j in range(i+1, nmols):
                pdtid = indices[j]
                pdtformula = formulas[j]
                pdt_elem_count = aseformula(pdtformula).count()
                if react_elem_count == pdt_elem_count:
                    reaction_ids.append([reactid, pdtid])
            if (i+1) % 1000 == 0:
                logging.info("Number of reactants iterated: %d" % (i+1))
                logging.debug("Number of rows to be writen to outfile: %s is %d" % (output, len(reaction_ids)))
                writer.writerows(reaction_ids)
                reaction_ids = []
        # Write if there is any reminder after the loop
        if reaction_ids:
            writer.writerows(reaction_ids)
    return


def get_arguments():
    parser = argparse.ArgumentParser("File to create bond list from xyz file.")
    parser.add_argument(
        "-mol_data", "--mol_data",
        type=str,
        required=True,
        help="csv file containing all QM9 molecules data"
    )
    parser.add_argument(
        "-output", "--output",
        type=str,
        required=False,
        default="reactions.csv",
        help="Output csv file for the reaction indices."
    )
    return parser.parse_args()


if __name__ == "__main__":
    st = time.time()
    args = get_arguments()
    molecule_csv_file = args.molecule_data
    output_csv_file = args.output
    #
    # End of parsing arguments
    #
    find_reactions(molecule_csv_file, output=output_csv_file)
    logging.info("Finding reaction completed.")
    logging.info("Job finished.")
    et = time.time()
    logging.info("Total time taken: %s" % round((et-st)/3600.0, 2))
