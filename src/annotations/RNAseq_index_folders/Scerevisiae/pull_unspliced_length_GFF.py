import argparse
from datetime import datetime
import sys

def calculate_unspliced_transcript_lengths(gff_file, output_tsv, command_line):
    cds_by_gene = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                feature_type = columns[2]
                start = int(columns[3])
                end = int(columns[4])
                attributes = columns[8]

                if feature_type == "CDS":
                    gene_ids = [x.split("=")[1] for x in attributes.split(";") if x.startswith("gene_id")]
                    if gene_ids:
                        gene_id = gene_ids[0]
                    else:
                        continue
                    if gene_id not in cds_by_gene:
                        cds_by_gene[gene_id] = []
                    cds_by_gene[gene_id].append((start, end))

    with open(output_tsv, 'w') as out:
        # Write the header with date, time, and command used
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        out.write(f"# Date and Time: {current_time}\n")
        out.write(f"# Command: {command_line}\n")
        out.write("Gene_ID\tunspliced.length\n")
        
        for gene_id, cds_coords in cds_by_gene.items():
            cds_coords = sorted(cds_coords, key=lambda x: x[0])
            first_cds_start = cds_coords[0][0]
            last_cds_end = cds_coords[-1][1]
            length = last_cds_end - first_cds_start + 1
            out.write(f"{gene_id}\t{length}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate unspliced transcript lengths from a GFF file.")
    parser.add_argument("gff_file", help="Path to the input GFF file.")
    parser.add_argument("output_tsv", help="Path to the output TSV file containing the lengths.")
    
    args = parser.parse_args()
    command_line = ' '.join(['python'] + [arg for arg in sys.argv])
    calculate_unspliced_transcript_lengths(args.gff_file, args.output_tsv, command_line)

if __name__ == "__main__":
    main()
