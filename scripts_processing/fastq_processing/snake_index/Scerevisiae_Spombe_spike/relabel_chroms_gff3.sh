#!/bin/bash

while read line; do
    if [[ -z "$line" ]] || [[ "$line" == \#* ]]; then
        # Skip empty lines and lines starting with #
        echo "$line"
    else
        # Prepend "pombe" to non-empty, non-comment lines
        echo "pombe$line"
    fi
done < "$1"
