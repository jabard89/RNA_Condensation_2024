library(tidyverse)
library(conflicted)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

github_dir <- here()
df_raw_samplesheet <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv"))

df_samples_byLysate <- df_raw_samplesheet %>%
  select("Lysate_sample","Strain","Treatment_group","Treatment","Treated","Temperature","Time","Biorep") %>%
  unique %>%
  write_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))

