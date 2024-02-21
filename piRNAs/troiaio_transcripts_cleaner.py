def process_fasta(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_header = None
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                current_header = line.split()[0]
                if 'type=' in line:
                    type_info = line.split('type=')[1].split(';')[0]
                    current_header += '_' + type_info
                f_out.write(current_header + '\n')
            else:
                f_out.write(line + '\n')

if __name__ == '__main__':
    input_file = '/Users/ascarpa/Downloads/dsim-all-miRNA-r2.02.fasta'
    output_file = '/Users/ascarpa/Downloads/dsim-all-miRNA-clean.fasta'
    process_fasta(input_file, output_file)

