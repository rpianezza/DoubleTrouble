import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Rename sequences in fasta file.")
parser.add_argument("fasta", help="Path to fasta file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

'''
Takes as input a fasta file and renames the TE insertions in a subsequential order by taxa.
es. Dmel-1, Dmel-2, Dmel-3 (if Dmel is the basename of the file)
'''

def rename(path, out):
    basename = os.path.splitext(os.path.basename(path))[0]
    genre = basename.split('.')[0]
    species = basename.split('.')[1][0:3]
    if len(basename.split('.'))>2:
        strain_start = basename.split('.')[2][0]
        strain_end = basename.split('.')[2][-1]
        taxa = str(genre)+str(species)+str(strain_start)+str(strain_end)
    else:
        taxa = str(genre)+"."+str(species)
    n = 1 
    with open(path, 'r') as input_file, open(out, 'w') as output_file:
        for line in input_file:
            if line.startswith(">"):
                name = ">"+str(taxa)+"-"+str(n)
                n+=1
                output_file.write(name+"\n")
            else:
                output_file.write(line)
            
rename(args.fasta, args.output)
