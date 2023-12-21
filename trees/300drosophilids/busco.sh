#!/bin/bash

# Set input and output directories
input_directory="/Volumes/Storage/data/assemblies-droso-petrov/"
output_directory=""
destination_directory="/Volumes/EXT-RICCARDO/DoubleTrouble/species-tree/"

# Iterate over all *fa files in the input directory
for fa_file in "${input_directory}"*.fasta; do
    # Get the assembly basename without extension
    assembly_basename=$(basename "${fa_file}" .fasta)

    # Create a subfolder in the output directory
    assembly_output_directory="${output_directory}${assembly_basename}/"
    mkdir -p "${assembly_output_directory}"

    # Run the busco command
    busco -i "${fa_file}" -l diptera_odb10 -o "${assembly_output_directory}" -m genome -c 20 -f

    # Move the subfolder to the destination directory
    mv "${assembly_output_directory}" "${destination_directory}"
done
