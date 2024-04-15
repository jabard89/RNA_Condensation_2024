# Load kallisto files
library(tidyverse)
library(vroom)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
df_samples <- read_tsv("RNAseq_samplesheet_with_kallisto_230510.tsv")
github_dir <- c("../../")
gene_labels <- read_tsv(paste0(github_dir,"/src/annotations/230508_labeled_genes_scer.tsv"))
dataset_folders <- list.files("counts")
df_counts <- dataset_folders %>% map(function(dataset){
  file_names <- list.files(file.path("counts",dataset,"counts/counts"),full.names=F)
  file_paths <- list.files(file.path("counts",dataset,"counts/counts"),full.names=T)
  file_paths %>% map2(file_names,function(file,name){
    read_tsv(file,show_col_types = F,comment="#") %>%
      mutate("RNA_sample"=str_replace(name,"_counts.tsv",""),
             "Dataset"=dataset)
  }) %>% bind_rows
}) %>% bind_rows %>%
  left_join(df_samples,by=c("RNA_sample","Dataset")) %>%
  select(-c(RNA_sample,Kallisto_file))

sedseq_mixing_ratios <- read_tsv("sedseq_mixing_ratios_230510.tsv")

ratios_join <- sedseq_mixing_ratios %>% filter(Term=="x50",Variable!="phi",Variable!="lp__") %>%
  mutate(Fraction=str_extract(Variable,"(?<=mixing_).*")) %>%
  rename("Mixing_ratio"="Value") %>%
  select(Dataset,Lysate_sample,Fraction,Mixing_ratio)

df_pSup <- df_counts %>%
  filter(Dataset!="snake_220328",
         !(Dataset=="snake_220725"&Lysate_sample=="L16"),
         !(Dataset=="snake_230428"&Lysate_sample=="L88")) %>%
  filter(Experiment=="SedSeq") %>%
  mutate(Total.counts = intron+unspliced+spliced+other) %>%
  pivot_longer(cols=c(intron,unspliced,spliced,other,Total.counts),names_to="Count_type",values_to="Counts") %>%
  filter(Counts>=20) %>%
  left_join(ratios_join,by=c("Lysate_sample","Fraction","Dataset")) %>%
  mutate(Mixing_ratio=if_else(Fraction=="Total",1,Mixing_ratio)) %>%
  mutate(Counts.adj = Counts*Mixing_ratio) %>%
  select(-c(Counts,Mixing_ratio,TPM)) %>%
  pivot_wider(names_from=Fraction,values_from=Counts.adj) %>%
    mutate(pSup = Sup/(Sup+Pellet))
df_Zsup_kallisto <- read_tsv("sedseq_comp_duplicates_230510.tsv.gz")

gene_labels_spliced <- gene_labels %>%
  filter(LengthTxEst!=LengthTxUnsplicedEst) %>%
  mutate(Spliced=FALSE,
         LengthTxEst=LengthTxUnsplicedEst)
gene_labels_all <- gene_labels %>%
  mutate(Spliced=TRUE) %>%
  bind_rows(gene_labels_spliced)

odds <- function(p) {
  p/(1-p)
}
logodds <- function(x) {
  # log odds, a shortcut
  res <- log(odds(x))
  res[!is.finite(res)] <- NA
  res
}
calc_mean_sd_roll <- function(df, col_name) {

  total_counts <- df %>%
    filter(Count_type == "Total.counts", Spliced == TRUE) %>%
    ungroup %>%
    select(LengthTxEst, all_of(col_name)) %>%
    arrange(LengthTxEst)

  means <- numeric(nrow(df))
  sds <- numeric(nrow(df))

  for (i in 1:nrow(df)) {
    current_length <- df$LengthTxEst[i]

    # Separate data into larger and smaller LengthTxEst values
    larger <- head(total_counts %>% filter(LengthTxEst >= current_length), 50)
    smaller <- tail(total_counts %>% filter(LengthTxEst < current_length), 50)
    subset_rows <- rbind(smaller, larger)

    # Calculate mean and sd for the given column
    means[i] <- mean(subset_rows[[col_name]], na.rm = TRUE)
    sds[i] <- sd(subset_rows[[col_name]], na.rm = TRUE)
  }

  result <- df %>%
    mutate(!!paste0("mean.", col_name) := means,
           !!paste0("sd.", col_name) := sds)

  return(result)
}

df_Zsup <- df_pSup %>%
  mutate(Spliced=if_else(Count_type%in%c("other","Total.counts","spliced"),TRUE,FALSE)) %>%
  group_by(Lysate_sample,ORF) %>%
  left_join(gene_labels_all,by=c("ORF","Spliced")) %>%
  filter(classification=="Verified") %>%
  filter(!is.na(pSup)) %>%
  mutate(lo.pSup=logodds(pSup)) %>%
  group_by(Dataset,Lysate_sample) %>%
  nest %>%
  mutate(data=map(data,function(df){
    print(df)
    # calculate Zsup
    # Sort by length so bins follow length
    df_Zsup <- calc_mean_sd_roll(df,"lo.pSup") %>%
      mutate(Zsup = (lo.pSup-mean.lo.pSup)/sd.lo.pSup)
    return(df_Zsup)
  })) %>%
  unnest(data)
write_tsv(df_Zsup,"all_counts_pSup_Zsup.tsv.gz")


