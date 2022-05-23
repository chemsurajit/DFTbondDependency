# DFTbondDependency
1) Repository for bond dependency paper

Dependency:
1) RDKit
2) xyz2mol
3) pandas
4) csv
5) Numpy

To install the xyz2mol package, consult the webpage:https://github.com/jensengroup/xyz2mol.git

All the other packages are part of standard Python library and can be installed with either Conda or PIP.


Description of the Data.
The data is part of the SURE project and can be downloaded from: <link>

Description of the codes.
This repository contains three python scripts to be run serially.
The scripts and their purpose is explained below:
1) get_bond_data.py
This script will download the data from the website <link> as zip file. Then, it will extract the zip files into the directory csvs.
Three type of files will be downloaded.
a) The log files from SCM calculation for all the molecules in energy refined QM9 database (link). The extracted data will be saved into the atoms and molecules subdirectories and the corresponding csv files atoms.csv and molecules.csv file will be stored.
b) The <qm9_molecules_bonds_energies.csv> file which contains the bonds and energies of the QM9 molecules.
c) The <reactions.csv> file which contains the indices of the reactant and products. The indices indicates the index column of the molecule dataset.

2) make_qm9_mol_bonds_ens.py
This script will create the csv files described above from the log files inside "TZP" directory.
3) make_reactions_bonds_ens.py
This script will create the csv file for the bond changes and reaction energies from the qm9_mols.csv, and reactions.csv files.
It additionally neads a csv file containing the QM9 molecular indices and the G4MP2 energy values (Which is provided as a supporting info in the paper: <link>)

4) do_lr_regression.py
This script will perform the linear regression described in paper <link> using the data and provide the coefficients in a csv file.

5) calculate_reaction_energy.py
This script will calculate the reaction energy from the LR coefficients found by the above scrpt.
It will correct the energy and provide the reaction energy in PBE/B3LYP-D/M06-2X and in G4MP2 level of theory.


How to run:
