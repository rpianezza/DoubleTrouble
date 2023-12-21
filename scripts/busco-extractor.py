# Takes as input the folder where the BUSCO subfolders are stored (output of "busco.sh").
# Outputs a fasta file for each BUSCO entry shared by all the species and a list of the shared BUSCO entries 

import os
import argparse

def read_entries(file_path):
    entries = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if line.startswith("#"):
                pass
            elif columns[1] == 'Complete':
                entries.add(columns[0])
    return entries

def find_common_entries(root_folder):
    common_entries = set()

    for subfolder in os.scandir(root_folder):
        if subfolder.is_dir():
            subsubfolder_path = os.path.join(subfolder.path, 'run_diptera_odb10')
            if os.path.exists(subsubfolder_path):
                current_entries = read_entries(os.path.join(subsubfolder_path, 'full_table.tsv'))
                if not common_entries:
                    common_entries.update(current_entries)
                else:
                    common_entries.intersection_update(current_entries)

    return list(common_entries)

def extract_sequences(gene, root_folder, output):
    for subfolder in os.scandir(root_folder):
        if subfolder.is_dir():
            subsubfolder_path = os.path.join(subfolder.path, 'run_diptera_odb10')
            if os.path.exists(subsubfolder_path):
                fasta_file_path = os.path.join(subsubfolder_path, 'busco_sequences', 'single_copy_busco_sequences', f'{gene}.faa')
                if os.path.exists(fasta_file_path):
                    species_name = subfolder.name
                    print(species_name)
                    with open(fasta_file_path, 'r') as fasta_file, open(output, 'a') as output_file:
                        for record in fasta_file:
                            if record.startswith('>'):
                                output_file.write('>'+species_name+'\n')
                            else:
                                output_file.write(record)

def main():
    parser = argparse.ArgumentParser(description='Extract and find common entries from TSV files.')
    parser.add_argument('input_folder', help='Path to the input folder containing subfolders with TSV files')
    parser.add_argument('output_folder', help='Path to the output file where the common entries will be saved')

    args = parser.parse_args()

    result = find_common_entries(args.input_folder)

    with open(args.output_folder+'/'+"shared-busco.txt", 'w') as output_file:
        for entry in result:
            output_file.write(entry + '\n')

    for gene in result:
        print(gene)
        extract_sequences(gene, args.input_folder, args.output_folder+'/'+gene+'.fasta')
        print("\n")

if __name__ == "__main__":
    main()