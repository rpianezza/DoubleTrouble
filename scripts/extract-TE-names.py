import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Rename sequences in phylip file.")
parser.add_argument("fasta", help="Path to MSA.fasta file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

'''
Takes as input a fasta library, output a txt file with every sequence name in each row.
Remove parts of the name after white space.
'''
def rename(fasta, out):
    with open(fasta, 'r') as fasta:
        with open(out, 'w') as output_file:
            for line in fasta:
                if line.startswith(">"):
                    name = line.strip().split(' ')[0]
                    new_name = str(name[1:]) + "\n" 
                    output_file.write(new_name)
                else:
                    continue

rename(args.fasta, args.output)
