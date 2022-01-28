from ase.formula import Formula as aseformula
import argparse
import pandas as pd
import csv
import time


def find_reactions(mol_csv_data, output="", st=0):
    """function to get reactions."""
    indices = mol_csv_data["index"].tolist()
    smiles = mol_csv_data["smiles"].tolist()
    formulas = mol_csv_data["chemformula"].tolist()
    nmols = len(indices)
    print("Reaction search will start now.")
    with open(output, 'w') as fcsv:
        writer = csv.writer(fcsv, delimiter=",", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerow(['reactindex', 'pdtindex'])
        for i in range(nmols):
            reactid = indices[i]
            reactformula = formulas[i]
            react_elem_count = aseformula(reactformula).count()
            for j in range(i+1, nmols):
                pdtid = indices[j]
                pdtformula = formulas[j]
                pdt_elem_count = aseformula(pdtformula).count()
                if react_elem_count == pdt_elem_count:
                    writer.writerow([reactid, pdtid])
            if (i+1) % 1000 == 0:
                print("Number of outerloop iterated: ", (i+1))
                print("Time: %f hrs" % ((time.time()-st)/3600))
    return


def main():
    #
    # parsing argument
    #
    st = time.time()
    parser = argparse.ArgumentParser("File to create bond list from xyz file.")
    parser.add_argument('-m', '--molecule_data', help='Location of the qm9 data extracted.',
                        required=True)
    parser.add_argument('-o', '--output', help='Output csv file for the reaction',
                        required=False, default='reactions.csv')
    args = parser.parse_args()
    molecule_csv_file = args.molecule_data
    output_csv_file = args.output
    #
    # End of parsing arguments
    #
    mol_csv_data = pd.read_csv(molecule_csv_file, usecols=["index", "smiles", "chemformula"])
    find_reactions(mol_csv_data, output=output_csv_file, st=st)
    return


if __name__ == "__main__":
    main()
