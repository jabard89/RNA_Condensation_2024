#!/bin/bash

# Check if a filename is provided as an argument
if [[ $# -eq 0 ]]; then
    echo "Usage: bash calculate_gc_percentage.sh <filename.fasta>"
    exit 1
fi

filename=$1

# Check if the file exists
if [[ ! -f $filename ]]; then
    echo "File not found: $filename"
    exit 1
fi

# Function to calculate GC percentage
calculate_gc_percentage() {
    local sequence=$1
    local length=${#sequence}
    local gc_count=$(echo -n "$sequence" | tr -cd 'GCgc' | wc -m)
    local gc_percentage=$(echo "scale=2; $gc_count * 100 / $length" | bc)
    echo $gc_percentage
}

# Read the contents of the file
sequence=""
while IFS= read -r line; do
    if [[ $line =~ ^\>.* ]]; then # we reached the next header
        if [[ ! -z $sequence ]]; then
            # Calculate the GC percentage for the sequence
            gc_percentage=$(calculate_gc_percentage "$sequence")

            # Output the sequence name and GC percentage as TSV
            echo -e "$sequence_name\t$gc_percentage"
        fi
        header_line=${BASH_REMATCH[0]}
        sequence_name=$(echo "$header_line" | awk '{print $1}' | tr -d '>')
        sequence=""
    else
        # Concatenate the sequence lines without spaces
        sequence+=$(echo "$line" | tr -cd 'ATGCUatgcu')
    fi
done < "$filename"

# Calculate GC percentage for the last sequence
gc_percentage=$(calculate_gc_percentage "$sequence")

# Output the sequence name and GC percentage as TSV
echo -e "$sequence_name\t$gc_percentage"
