#!/bin/bash
#SBATCH --job-name={SAMPLE}
#SBATCH --output=sbatch/map_{SAMPLE}_%j.out
#SBATCH --error=sbatch/map_{SAMPLE}_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-dadrummond
#SBATCH --mem=32G

module load python/anaconda-2022.05
source activate /home/jbard/beagle3-dadrummond/jbard/envs/py310_snake_star
cd {BASEDIR}
snakemake --snakefile {SNAKEFILE} --config sample={SAMPLE} \
fastq1="{FASTQ1}" \
STAR_index_dir={STAR} \
gtf={GTF} \
kallisto_index={KALLISTO} \
kallisto_direction="{KALLISTO_DIRECTION}" \
--cores all -p --use-conda
