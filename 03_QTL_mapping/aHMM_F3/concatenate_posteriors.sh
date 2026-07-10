#!/bin/bash

# Define the directories
dir1="chr1"
dir2="chr2"
dir3="chr3"
dir4="chr4"
dir5="chr5"
dir6="chr6"
dir7="chr7"
dir8="chr8"
dir9="chr9"

# Output directory
output_dir="concatenated_files_unthinned"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through files in the first directory
for file in "$dir1"/*; do
    filename=$(basename "$file")
    # Remove header and concatenate files
    tail -n +2 "$dir1/$filename"  > "$output_dir/$filename"
    tail -n +2 "$dir2/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir3/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir4/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir5/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir6/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir7/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir8/$filename"  >> "$output_dir/$filename"
    tail -n +2 "$dir9/$filename"  >> "$output_dir/$filename"
done

echo "Files concatenated successfully!"