#!/bin/bash
fold_dir="/home/jabard89/Dropbox/code_JB/repos/RNAfold"
pyscript_dir="/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools"
python ${pyscript_dir}/scripts/annotation/extract_UTR_from_GFF.py \
--add_nts 20 \
${pyscript_dir}/UTRs/Saccharomyces_cerevisiae.R64-1-1.109_Pelechano2013.gff3 \
${pyscript_dir}/UTRs/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
${fold_dir}/input/Scerevisiae_UTR5plus20.csv

grep -vwE "(ORF|^#)" ${fold_dir}/input/Scerevisiae_UTR5plus20.csv > \
${fold_dir}/input/Scerevisiae_UTR5plus20_clean.csv

output_dir=${fold_dir}/output
while IFS="," read ORF UTR5 UTR5_plus20 UTR3; do
	if [ -n $UTR5_plus20 ]; then
		mkdir -p $output_dir/ORFs/"${ORF}"
		printf ">${ORF}\n${UTR5_plus20}" > "$output_dir/ORFs/${ORF}/${ORF}_UTR5plus20.fasta"
	fi
done < ../input/Scerevisiae_UTR5plus20_clean.csv


for dir in $output_dir/ORFs/*; do
	orf=$(basename $dir)
	RNAfold -p0 --noPS --noDP \
		--infile="$dir/${orf}_UTR5plus20.fasta" \
		> "$dir/${orf}_UTR5plus20.fold"
	echo "Finished $orf"
done


output_file=$output_dir/ViennaFold_UTR5plus20.tsv
printf "%b\n" "ORF\tEnergy" > ${output_file}
for dir in $output_dir/ORFs/*; do
	orf=$(basename $dir)
	energy=$(grep -P -o -e '[-0-9.]*(?= kcal)' $dir/${orf}_UTR5plus20.fold)
	printf "%b\n" "${orf}\t${energy}" >> ${output_file}
done