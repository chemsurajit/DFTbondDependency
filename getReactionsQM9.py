from ase.formula import Formula as aseformula
import argparse
import pandas as pd
import csv



def find_reactions(mol_csv_data, output=""):
    """function to get reactions."""
    indices = mol_csv_data["index"].tolist()
    smiles = mol_csv_data["smiles"].tolist()
    formulas = mol_csv_data["chemformula"].tolist()
    nmols = len(indices)
    pd_react = pd.DataFrame(columns=['reactindex', 'reactsmi', 'pdtindex', 'pdtsmi'])
    with open(output,'w') as fcsv:
        writer = csv.writer(fcsv, delimiter=",", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerow(['reactindex','reactsmi','pdtindex','pdtsmi'])
        for i in range(nmols):
            reactid = indices[i]
            reactsmi = smiles[i]
            reactformula = formulas[i]
            react_elem_count = aseformula(reactformula).count()
            for j in range(i+1, nmols):
                pdtid = indices[j]
                pdtsmi = smiles[j]
                pdtformula = formulas[j]
                pdt_elem_count = aseformula(pdtformula).count()
                if react_elem_count == pdt_elem_count:
                    writer.writerow([reactid,reactsmi,pdtid,pdtsmi])
            if (i+1) % 10000 == 0:
                print("Number of outerloop iterated: ", (i+1))
    return


def main():
    #
    # parsing argument
    #
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
    mol_csv_data = pd.read_csv(molecule_csv_file, usecols=["index", "smiles", "chemformula"], index=False, nrows=100)
    find_reactions(mol_csv_data, output=output_csv_file)

if __name__ == "__main__":
    main()
