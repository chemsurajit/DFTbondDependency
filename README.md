# DFTbondDependency
Repository for bond dependency paper

Dependency:
1) RDKit
2) xyz2mol
3) pandas
4) csv
5) Numpy
6) requests
7) statsmodel
8) sklearn

To install the xyz2mol package, visit: https://github.com/jensengroup/xyz2mol.git

All the other packages are part of standard Python library and can be installed with either Conda, PIP, etc.

If this repository is used, please cite us. Citation format can be downloaded by clicking the link: 
"Cite this repository" in the right side panel.

## Description of the Data.
The data is not publicly available now. It will be publicly available upon acceptance of our paper in a peer-reviewed journal. The preprint can be downloaded from: <link>

## Description of the codes.
The repository contains python scripts for the calculations described in the paper: <link after submission>
1) get_data.py: This script will download data from the DTU Data website. The public link to the website
will be available upon acceptance of the paper. This file will download either all the files from the 
database (if -all/--all option is given) Or it will only download the xyzfiles and the log files of the
energy calculations

2) make_molecule_bond_en_csv.py: This file makes a csv file containing energy values, list of bonds, SMILES
string, chemical formula, etc from the xyzfiles and logfiles downloaded by using the get_data.py script.

3) make_reaction_ids.py: This script will create a csv file with only two columns: 'reactantindex','pdtindex'.
The indices are the index of the molecules in the csv file made by the make_molecule_bond_en_csv.py script.

4) process_reaction_conversion_jobs.sh: This is a bash script to create the final csv file containing
all the information related to the reactions. It takes the csv file containing molecular data (created by the
script make_molecule_bond_en_csv.py), the indices of the reactants and products in form of a csv file
(created by using the script make_reaction_ids.py), The G4MP2 energies of the molecules as csv file (with
index and energy), the path of the python script make_reactions_parallel.py, number of Nodes to be used,
and number of processors per each nodes.
It first split the csv file containing indices of the reactants and products according to the number of Nodes
and saves those in a json file with names Node_n.json with n from {1,2,...n} if n number of nodes are used.

5) make_reactions_parallel.py: This file takes csv file containing indices for the reactions, csv file containing
all the data of the molecules, csv file containing G4MP2 energy, number of processors, json file
containing the indices of the csv file with "reactantindex","pdtindex".

6) submit.sh: An example submit script to run make_reactions_parallel.py in a single node with multiple processors.
It is called from the script process_reaction_conversion_jobs.sh. It is written for the slurm scheduler.

7) detect_correlation.py: This script is for detecting correlation between the variables (bonds). It takes
the directory location of the csv files containing all the reaction data (created by the process_reaction_conversion_jobs.sh)
script. By default, it randomly chose 10% of the total data to detect correlation.

8) do_linear_regression.py: This script performs the linear regression between the bond change and the DFT error to
reaction energies. It takes as argument the directory location for the reaction data file, and the name of the DFT
functionals.

9) correct_reaction_energy.py: This script calculates the reaction energy and the correction to it for a given DFT functional.
It takes as input the log files of reactants (with the option -r), products (with the option -p), and name of the
DFT functional and prints out the reaction energy for the DFT functional, the correction, and the corrected reaction energy.

10) CITATION.cff: This file is to provide citation data for this repository in bibtex or APA format.

## How to run:
The help message for each of the files (except submit.sh) can be obtained by running the corresponding
script with -h.


The steps described in the paper can be followed by running the below scripts in the following sequence: 
1. get_data.py
2. make_molecule_bond_en_csv.py
3. make_reaction_ids.py
4. process_reaction_conversion_jobs.sh, make_reactions_parallel.py, submit.sh
5. detect_correlation.py
6. do_linear_regression.py
