import sys
import argparse
import csv
import pandas as pd
import time

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



def build_csv_for_reactions(start_line_no, nline, reaction_csv_file, molecule_csv_file, outputcsv):
    start_row = start_line_no
    end_row = start_row + nline
    pd_reactions = pd.read_csv(reaction_csv_file, header=0, skiprows=range(1,start_row), nrows=nline )
    #sliced_rean_pd = pd_reactions.loc[start_row:end_row]
    pd_molecule = pd.read_csv(molecule_csv_file)
    bonds_energy_funcs_list = bonds_list + dft_functionals
    mode = 'w'
    rindex_for_csv = start_row
    #whole_df = pd.DataFrame()
    for index, row in pd_reactions.iterrows():
        loop_st = time.time()
        reactant_index = row.reactindex
        pdt_index = row.pdtindex
        reactant_row = pd_molecule.loc[pd_molecule['index'] == reactant_index]
        pdt_row = pd_molecule.loc[pd_molecule['index'] == pdt_index]
        reactant_smi = reactant_row["smiles"].values[0]
        pdt_smi = pdt_row["smiles"].values[0]
        rtn_prop_diff = pdt_row[bonds_energy_funcs_list] - reactant_row[bonds_energy_funcs_list].values
        rtn_prop_diff["recursive_id"], rtn_prop_diff["reactionindex"], rtn_prop_diff["reactantindex"], rtn_prop_diff["pdtindex"], rtn_prop_diff["react_smi"], rtn_prop_diff["pdt_smi"] = \
            [index, rindex_for_csv, reactant_index, pdt_index, reactant_smi, pdt_smi]
        #whole_df = whole_df.append(rtn_prop_diff)
        if mode == 'w':
            rtn_prop_diff.to_csv(outputcsv, mode=mode, sep=",", index=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
        if mode == 'a':
            rtn_prop_diff.to_csv(outputcsv, mode=mode, sep=",", index=False, header=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
        rindex_for_csv += 1
        mode = 'a'
        print("Loop timing: ", time.time()-loop_st)
        if (index%10000) == 0:
            print("Index: ", index)
    #whole_df.to_csv(outputcsv, mode='w', sep=",", index=False, quoting=csv.QUOTE_MINIMAL, na_rep='nan')
    return


def main():
    #
    # parsing arguments
    #
    parser = argparse.ArgumentParser("File to generate reactions with the energy, bonds etc.")
    parser.add_argument('-r', '--reaction_csv', help="CSV file location (reactions.csv)", required=False)
    parser.add_argument('-m', '--molecule', help="CSV file containing molecular bond and energy data", required=False)
    parser.add_argument('-o', '--output', help="File name for the output csv file.", required=False,
                        default="core.csv")
    parser.add_argument('-sl', '--start_line', help="Line number of the starting for reactions csv file.", required=False,
                        default=0)
    parser.add_argument('-nl', '--nline', help="How many lines to be read from the reactions csv file.", required=False,
                        default=1000)
    args = parser.parse_args()
    reaction_csv_file = args.reaction_csv
    #reaction_csv_file = "../reactions.csv"
    molecule_csv_file = args.molecule
    #molecule_csv_file = "../qm9_bonds_energies.csv"
    output_csv_file = args.output
    #output_csv_file = "core.csv"
    start_line_no = int(args.start_line)
    #start_line_no = int(sys.argv[1])
    nlines = int(args.nline)
    #nlines = int(sys.argv[2])
    #
    # parsing done
    #
    print("start line no: ", start_line_no)
    build_csv_for_reactions(start_line_no, nlines, reaction_csv_file, molecule_csv_file, output_csv_file)
    print("Job finished.")
    print("-------------")
    return 


if __name__ == "__main__":
    main()
