# Loop_extractor

It is a small tool to extract loops area between domains of various protein as linkers for fusion protein design etc. (Python 3.7)

Usage:

1. Change the line 34 to the applied protein dataset path

3. Run the program, it will read the uniprot id of the input protein dataset, and extract related information automatically

# Loop_selector

For any input .pdb file, the program will extract loop area of this protein utilizing pyrosetta loop selector funciton, and calculate the related information using GOR4 method, and storge in an excel file

This can help to select loop area of a natural protein as protein linker design

Usage:

change the runpath in pdb_scan.py to the file with all protein needed to be extract
Then the program will run the loop_selector.py to extract loop area of all the protien and store all the information collected in an excel file
