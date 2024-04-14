# Calculate Mixing Ratios
# Jared Bard
# 230504 updated 240414
library(tidyverse)
library(rstan)
library(conflicted)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(rstan::extract)
conflicts_prefer(stats::lag)
conflicts_prefer(purrr::compose)
conflicts_prefer(cowplot::align_plots)

# Specify a seed for random number generator so output is reproducible
myseed = 40
#run in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

github_dir <- here()
source(file.path(github_dir,"scripts_processing/ewj_mixing_functions.R"))

gene_labels <- read_tsv(file.path(github_dir,"src/annotations/labeled_genes_scer.tsv")) 

df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))) %>%
  mutate(Dataset = paste0("snake_",Sequencing_date))

read_kallisto <- function(file) {
  kallisto_coltypes <- list(ORF="c",length="i",eff_length="n",est_counts="n",TPM="n")
  vroom::vroom(file,col_names=c("ORF","length","eff_length","est_counts","TPM"),
        col_types=kallisto_coltypes,
        skip=1) %>%
    mutate("Kallisto_file"=file)
}

prepare_sedseq_data_for_fitting <- function(df) {
  df %>%
    left_join(df_samples,by="Kallisto_file") %>%
    left_join(gene_labels,by="ORF") %>%
    dplyr::filter(classification=="Verified") %>%
    filter(Fraction%in%c("Total","Sup","Pellet")) %>%
    group_by(Dataset,Lysate_sample) %>%
    dplyr::filter(est_counts>quantile(est_counts,0.05)) %>%
    mutate(round_counts = round(est_counts)) %>%
    ungroup %>%
    select(ORF,Dataset,Lysate_sample,round_counts,Fraction) %>%
    pivot_wider(names_from="Fraction",values_from="round_counts") %>%
    drop_na()
}

df_sedseq_example <- df_samples %>%
  filter(Lysate_sample=="L7") %>%
  pull(Kallisto_file) %>%
  map(function(file) {read_kallisto(file)}) %>% bind_rows
sedseq_columns <- c("Sup","Pellet")
# compile stan model, test it with 10 iterations, save output

sedseq_stan_mixing_model <- initialize_stan_mixing_model(prepare_sedseq_data_for_fitting(df_sedseq_example),sedseq_columns)

# run stan mixing ratio inference, reusing same model
sedseq_mixing_ratios <- df_samples %>%
  filter(Experiment=="SedSeq") %>%
  group_by(Dataset,Lysate_sample) %>% nest %>%
  mutate(data=pmap(list("df"=data,"dataset"=Dataset,"Lysate_sample"=Lysate_sample),function(df,dataset,Lysate_sample) {
    file_out <- file.path(github_dir,"data_processed/mixing_ratios/",
                       paste0(dataset,"~",Lysate_sample,".tsv"))
    if (file.exists(file_out)){
      return(read_tsv(file_out,show_col_types=F))
    }
    else {
      print(df$Kallisto_file)
      df_raw <- df %>% pull(Kallisto_file) %>%
        map(function(file) {read_kallisto(file)}) %>% bind_rows
      
      mod <- fit_mixing_model(prepare_sedseq_data_for_fitting(df_raw),sedseq_columns,
                              stan_mixing=sedseq_stan_mixing_model,seed=myseed)
      ratios <- get_mixing_params(mod)
      write_tsv(ratios,file_out)
      return(ratios)
    }
  })
  ) %>% unnest(data)
write_tsv(sedseq_mixing_ratios,file.path(github_dir,"data_processed/sedseq_mixing_ratios.tsv"))
