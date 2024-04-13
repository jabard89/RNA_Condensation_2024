# Load the google spreadsheet
library(googlesheets4)
library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
#Read google sheets data into R
df_raw_samplesheet <- read_sheet('https://docs.google.com/spreadsheets/d/12GOWDh96JDPCGnYEeCU5uqhAm3fW3k_kOBHLDmLk104/edit?usp=sharing',
                                 sheet="Collated") %>%
  write_tsv("RNASeq_sample_sheet_230510.tsv")

df_samples_byLysate <- df_raw_samplesheet %>%
  select("Lysate_sample","Strain","Treatment_group","Treatment","Treated","Temperature","Time","Biorep") %>%
  unique %>%
  write_tsv("RNAseq_samples_byLysate_230510.tsv")

df_labels <- read_sheet('https://docs.google.com/spreadsheets/d/12GOWDh96JDPCGnYEeCU5uqhAm3fW3k_kOBHLDmLk104/edit?usp=sharing',
                                 sheet="Stress_labels") %>%
  mutate(Stress=factor(Stress,levels=c("none","Heat Shock","Azide","Azide pH4","Ethanol","Dextrose","DTT","KCl")),
         Stress_label=factor(Stress_label,levels=c("none","30C","42C","42CR","46C",
                                                   "mock","Azide 0.5%","Azide 0.8%",
                                                   "Ethanol 5%","Ethanol 7.5%","Ethanol 10%","Ethanol 15%",
                                                   "DTT 10mM","KCl 1M","withdrawal","Azide 0.5% pH4")))
write_tsv(df_labels,"stress_samples_230919.tsv")
