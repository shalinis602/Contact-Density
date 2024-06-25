import os
import warnings
from Bio import PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def read_coordinates(filename):
    """Read CA atom coordinates from a PDB file and return lists of x, y, z coordinates."""
    x_coords, y_coords, z_coords = [], [], []
    parser = PDB.PDBParser()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        structure = parser.get_structure('structure', filename)
    for atom in structure.get_atoms():
        if atom.get_id() == 'CA':
            x, y, z = atom.get_coord()
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)
    return x_coords, y_coords, z_coords, structure

def calculate_contacts(x_coords, y_coords, z_coords, distance_threshold=6.0):
    """Calculate contacts within the protein based on a distance threshold."""
    num_atoms = len(x_coords)
    contacts = []
    distance_threshold_squared = distance_threshold ** 2

    for i in range(num_atoms):
        for j in range(i + 2, num_atoms):  # Avoid self-contact and consecutive atoms
            squared_distance = (
                (x_coords[i] - x_coords[j]) ** 2 +
                (y_coords[i] - y_coords[j]) ** 2 +
                (z_coords[i] - z_coords[j]) ** 2
            )
            if squared_distance < distance_threshold_squared:
                contacts.append(1)
            else:
                contacts.append(0)
    return contacts

def calculate_common_contacts(contacts1, contacts2):
    """Calculate the number of common contacts between two contact lists."""
    common_contacts = sum(1 for c1, c2 in zip(contacts1, contacts2) if c1 and c2)
    return common_contacts

def align_structures(mobile_structure, target_structure):
    """Align the mobile structure to the target structure and return aligned coordinates."""
    super_imposer = PDB.Superimposer()
    mobile_atoms = [atom for atom in mobile_structure.get_atoms() if atom.get_id() == 'CA']
    target_atoms = [atom for atom in target_structure.get_atoms() if atom.get_id() == 'CA']
    
    super_imposer.set_atoms(target_atoms, mobile_atoms)
    super_imposer.apply(mobile_structure.get_atoms())
    
    x_coords, y_coords, z_coords = [], [], []
    for atom in mobile_atoms:
        x, y, z = atom.get_coord()
        x_coords.append(x)
        y_coords.append(y)
        z_coords.append(z)
    
    return x_coords, y_coords, z_coords

def main():
    # Prompt user for filenames
    conformation_filename = input("Enter the filename for the conformation structure (e.g., conformation.pdb): ")
    native_filename = input("Enter the filename for the native structure (e.g., native.pdb): ")

    # Read coordinates from the native conformation file
    x_coords_native, y_coords_native, z_coords_native, native_structure = read_coordinates(native_filename)

    # Read coordinates from the conformation file
    x_coords_conformation, y_coords_conformation, z_coords_conformation, conformation_structure = read_coordinates(conformation_filename)

    # Align the conformation structure to the native structure
    x_coords_conformation, y_coords_conformation, z_coords_conformation = align_structures(conformation_structure, native_structure)

    # Calculate contacts for the conformation file
    conformation_contacts = calculate_contacts(x_coords_conformation, y_coords_conformation, z_coords_conformation)
    conformation_contact_count = sum(conformation_contacts)

    # Calculate contacts for the native conformation file
    native_contacts = calculate_contacts(x_coords_native, y_coords_native, z_coords_native)
    native_contact_count = sum(native_contacts)

    # Calculate the number of common contacts
    common_contact_count = calculate_common_contacts(conformation_contacts, native_contacts)

    # Calculate contact density
    contact_density = (common_contact_count / native_contact_count) * 100 if native_contact_count > 0 else 0

    # Print results
    print("Contacts in given protein conformation =", conformation_contact_count)
    print("Contacts in native protein conformation =", native_contact_count)
    print("Common contacts =", common_contact_count)
    print("Contact density % =", contact_density)

if __name__ == "__main__":
    main()

