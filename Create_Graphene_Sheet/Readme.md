## Creating GNP sheets for PCFF
1. Start with GNP_lattice.in to create the GNP shape and size you want.
2. Visualize the data file in Ovito to verify that it's the shape and size you want.
3. Use the generate_graphitic_structure.py script with the data file from the previous step as input.
4. Use the convert_cnt_gnp___opls_2_iff.py script to convert the OPLS data file into a PCFF-IFF file with dummy atoms.
5. Visualize the final sheet in Ovito to verify it worked.
