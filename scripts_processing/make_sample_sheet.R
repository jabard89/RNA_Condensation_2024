library(tidyverse)
library(conflicted)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

github_dir <- here()
df_raw_samplesheet <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv"))
control_lysate <- c("L7","L11","L41","L44","L47","L48","L51","L54","L57","L63","L69","L75",
                    "L81","L83","L89","l90","L91","L92","L93","L94","L95",
                    "L96","L99","L102","L104","L106","L108","L110","L112","L114","L116","L122","L124")

df_samples_byLysate <- df_raw_samplesheet %>%
  select("Lysate_sample","Strain","Treatment_group","Treatment","Treated","Temperature","Time","Biorep") %>%
  unique %>%
  mutate(Control = ifelse(Lysate_sample %in% control_lysate,TRUE,FALSE)) %>%
  write_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))

