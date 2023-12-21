import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Rename sequences in FASTA file.")
parser.add_argument("fasta", help="Path to FASTA file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

'''
Takes as input a fasta library
'''
def rename(fasta, out):
    with open(fasta, 'r') as fasta_file:
        with open(out, 'w') as output_file:
            for line in fasta_file:
                if line.startswith(">"):
                    name = line.strip().split(' ')[0]
                    modified_name = ""
                    for letter in name:
                        if letter == '-':
                            letter = '_'
                        modified_name += letter
                    new_line = modified_name + "\n"
                    output_file.write(new_line)
                else:
                    output_file.write(line)

rename(args.fasta, args.output)