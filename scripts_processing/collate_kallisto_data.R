# Load kallisto files
library(tidyverse)
library(vroom)
library(conflicted)
library(here)

github_dir <- here()
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv")))

kallisto_coltypes <- list(ORF="c",length="i",eff_length="n",est_counts="n",TPM="n")

all_counts <- df_samples$RNA_sample %>% map(function(RNA_sample){
  file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))
  read_tsv(file,col_names=c("ORF","length","eff_length","est_counts","TPM"),
           col_types = kallisto_coltypes, skip = 1) %>% 
    mutate("Kallisto_file"=file,
           "RNA_sample"=RNA_sample)
}) %>% bind_rows

write_tsv(all_counts,file.path(github_dir,"data_raw","all_kallisto_counts.tsv.gz"))



