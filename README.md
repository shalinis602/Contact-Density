# Contact Density

This repository contains a Python script to calculate the contact density of a protein conformation with respect to its native conformation. It first aligns the given protein structure to the native structure and then calculates their contact density.

## Prerequisites

To run this script, you'll need to have the following installed:

- Python 3.6 or higher
- Biopython library

You can install Biopython using pip:

```sh
pip install biopython
```

## Usage

1. Clone the repository to your local machine:

```
git clone git@github.com:shalinis602/Contact-Density.git
cd Contact-Density
```

2. Place your protein conformation PDB files in the repository directory. Ensure you have two PDB files:

- One for the given protein conformation.
- One for the native protein conformation.

3. Run the script:

```
python contact_density.py
```

4. Follow the prompts to input the filenames of your protein conformation and native structure PDB files.

## Example

Here is an example of how to use the script:

```
Enter the filename for the conformation structure (e.g., conformation.pdb): markov_rf_4988.pdb
Enter the filename for the native structure (e.g., native.pdb): native.pdb
Contacts in given protein conformation = 100
Contacts in native protein conformation = 120
Common contacts = 80
Contact density % = 66.67
```

## How It Works

The script performs the following steps:

1. **Read Coordinates**: It reads the alpha carbon (CA) atom coordinates from the given PDB files.
2. **Align Structures**: It aligns the given protein conformation to the native structure using Biopython's `Superimposer`.
3. **Calculate Contacts**: It calculates the contacts within each structure based on a distance threshold (default is 6.0 Å).
4. **Calculate Common Contacts**: It calculates the number of common contacts between the two structures.
5. **Compute Contact Density**: It computes the contact density as the percentage of common contacts relative to the total contacts in the native structure.

## Contributing

If you find any issues or have suggestions for improvements or expanding the project, feel free to open an issue or submit a pull request.

## Acknowledgements

This script uses the Biopython library for PDB parsing and structure manipulation. Biopython is an open-source project available [here](https://biopython.org/).