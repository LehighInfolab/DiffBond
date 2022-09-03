# DiffBond

## Purpose:
DiffBond identifies and classifies three intermolecular bonds in protein complexes: ionic bonds, salt bridges, and hydrogen bonds. 


##Command-Line Options:

-i
Input PDB files to be compared if multiple. This option will default to finding intermolecular bonds at the interface between input PDBs.

-d
Resolution for distance checking. Increasing distance will lessen strictness in search, and vice versa. Default is 5 angstroms.

-m
Search mode option. Allows you to perform different searches instead of the default (all three of ionic bond, hydrogen bond, salt bridge search). Contact search = [], ionic bond search = i, hydrogen bond search =h, salt bridge search = s.

-o
Name of output file. Default is Contact_[PDB1]_[PDB2].pdb.


##Input:
The input can be two PDB files that form a complex and already in the correct spatial arrangement, or can be a single PDB file with target proteins indicated by different chain letters.

For hydrogen bond identification, input is a single PDB file with target proteins.


##Output:
The output consists of three lists of amino acid pairs: pairs that can form ionic bonds, pairs that can form hydrogen bonds, and pairs that can form salt bridges. Each list is an excel sheet in the format of [FirstPDB Chain AminoAcid SequenceNumbering] -> [Chain AminoAcid SequenceNumbering][Chain AminoAcid SequenceNumbering]..., where first column indicates target amino acid and all subsequent amino acids can form a bond with the first. Each pair will appear twice in the list, one in reverse order of the other.


##Options for visualization:
We suggest that users highlight relevant amino acids using their chain and sequence numbering from the output list. This can easily be done in PyMOL using the following command-line option: color red, resi [sequence number] and chain [chain letter].
