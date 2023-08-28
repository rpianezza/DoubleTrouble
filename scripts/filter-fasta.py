import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Filter FASTA file.")
parser.add_argument("fasta", help="Path to fasta file")
parser.add_argument("selected", help="Path to the file with the fasta sequences to select")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

def filter(fasta, selection, out):
    with open(fasta, 'r') as fasta_file, open(selection) as selected, open(out, 'w') as output_file:
        genes = []
        for line in selected:
            gene = line[0:-1]
            genes.append(gene)
        write_line = False
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                gene_name = line[1:]
                if gene_name in genes:
                    write_line = True
                    output_file.write(line + "\n")  # Write the current line to output
            elif write_line:
                output_file.write(line + "\n")  # Write the next line to output
                write_line = False

filter(args.fasta, args.selected, args.output)