#!/bin/bash

# Input file containing sample names
input_file="hyb_samps.list"

# Read each line from the file
while IFS= read -r sample_name; do
    # Create a file for each sample
    output_file="${sample_name}.txt"
    
    # Write the content to the file
    echo "parent1 19.467_VILL,19.468_VILL,19.470_VILL,19.475_VILL,19.476_VILL,19.582_VILL,19.583_VILL,19.584_VILL,19.585_VILL,19.587_VILL,19.588_VILL,19.600_VILL,19.601_VILL,19.602_VILL,19.603_VILL,19.604_VILL,19.605_VILL,19.728_VILL,19.729_VILL,21.234_VILL" > "$output_file"
    echo "parent2 19.495_ALLE,19.497_ALLE,19.500_ALLE,19.502_ALLE,19.504_ALLE,19.506_ALLE,19.508_ALLE,19.512_ALLE,19.513_ALLE,19.543_ALLE,19.544_ALLE,19.547_ALLE,19.552_ALLE,19.553_ALLE,19.555_ALLE,19.567_ALLE,19.596_ALLE,19.599_ALLE,19.606_ALLE,19.608_ALLE,19.609_ALLE,19.611_ALLE,19.613_ALLE,19.614_ALLE,19.720_ALLE,19.721_ALLE,19.722_ALLE,19.723_ALLE,19.724_ALLE,19.725_ALLE,19.730_ALLE" >> "$output_file"
    echo "admixed $sample_name" >> "$output_file"
done < "$input_file"

echo "Files generated successfully."
