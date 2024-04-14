#!/bin/bash
while read line; do
    if [[ "$line" == ">"* ]]; then
        # If the line starts with ">", it's a sequence ID line
        echo ">pombe${line#>}"  # add pombe to beginning of line (removing original >)
    else
        # Otherwise, it's a sequence line, so just echo it as-is
        echo "$line"
    fi
done < "$1"  # read input from the file specified as the first argument to the script
