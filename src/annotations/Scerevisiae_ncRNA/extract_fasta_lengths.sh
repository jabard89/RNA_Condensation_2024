#!/bin/bash

# Check that an input file was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 [input.fasta]"
    exit 1
fi

# Print the header line
echo -e "ORF\tlength"

# Loop through each sequence in the file
while read line; do
    # Check if the line is a header line
    if [[ $line == ">"* ]]; then
        # If it's a header line, print the previous sequence length (if any)
        if [ ! -z "$seq" ]; then
            echo -e "${name}\t${#seq}"
        fi
        # Save the new sequence name and reset the sequence string
        name=$(echo $line | cut -d " " -f 1 | sed 's/>//')
        seq=""
    else
        # If it's not a header line, append the sequence to the current string
        seq=${seq}${line}
    fi
done < "$1"

# Print the length of the last sequence
if [ ! -z "$seq" ]; then
    echo -e "${name}\t${#seq}"
fi
