# Jared Bard
# 230504
# use mixing ratios calculated by calculate_mixing_ratios to generate pSups and Zsups
library(tidyverse)
library(conflicted)
library(flextable)
library(zoo)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

github_dir <- here()
gene_labels <- read_tsv(file.path(github_dir,"src/annotations/labeled_genes_scer.tsv")) 
read_kallisto <- function(file) {
  kallisto_coltypes <- list(ORF="c",length="i",eff_length="n",est_counts="n",TPM="n")
  vroom::vroom(file,col_names=c("ORF","length","eff_length","est_counts","TPM"),
               col_types=kallisto_coltypes,
               skip=1) %>%
    mutate("Kallisto_file"=file)
}
df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))) %>%
  mutate(Dataset = paste0("snake_",Sequencing_date))
df_samples_byLysate <- read_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))
sedseq_mixing_ratios <- read_tsv(file.path(github_dir,"data_processed","sedseq_mixing_ratios.tsv"))
mixing_ratios <- bind_rows(sedseq_mixing_ratios)

#calculate the percent error in the fit to judge how tightly constrained the data is
mixing_ratios %>%
  filter(Term%in%c("mean","sd"),
         Variable!="phi",Variable!="lp__") %>%
  pivot_wider(names_from=Term,values_from=Value) %>%
  mutate(percent_error = sd/mean*100) %>%
  write_tsv(file.path(github_dir,"data_processed","mixing_ratio_summary.tsv"))

ratios_join <- mixing_ratios %>% filter(Term=="x50",Variable!="phi",Variable!="lp__") %>%
  mutate(Fraction=str_extract(Variable,"(?<=mixing_).*")) %>%
  rename("Mixing_ratio"="Value") %>%
  select(Dataset,Lysate_sample,Fraction,Mixing_ratio)

df_Zsup <- df_samples %>%
  left_join(ratios_join,by=c("Dataset","Lysate_sample","Fraction")) %>%
  filter(Experiment=="SedSeq") %>%
  select(Lysate_sample,Fraction,Kallisto_file,Mixing_ratio) %>%
  group_by(Lysate_sample) %>%
  nest %>%
  mutate(data=map(data,function(df){
    df_counts <- df$Kallisto_file %>% map(function(file) {read_kallisto(file)}) %>%
      bind_rows %>% left_join(df_samples,by="Kallisto_file") %>%
      left_join(df,by=c("Kallisto_file","Fraction"))
    df_tot <- df_counts %>% filter(Fraction=="Total") %>%
      mutate(Total.TPM = TPM,
             est_counts.Total=est_counts) %>%
      select(ORF,Total.TPM,est_counts.Total)
    # test <- df_counts %>% filter(Fraction!="Total") %>%
    #   select(ORF,Fraction,est_counts,Mixing_ratio) %>%
    #   mutate(counts_adj = est_counts*Mixing_ratio)
    df_pSup <- df_counts %>% filter(Fraction!="Total") %>%
      select(ORF,Fraction,est_counts,Mixing_ratio) %>%
      mutate(counts_adj = est_counts*Mixing_ratio) %>%
      pivot_wider(names_from=Fraction,values_from=c(est_counts,counts_adj,Mixing_ratio),names_sep=".") %>%
      mutate(mixSum=counts_adj.Sup+counts_adj.Pellet,
             pSup=counts_adj.Sup/mixSum,
             SP = counts_adj.Sup/counts_adj.Pellet) %>%
      left_join(df_tot,
                by=c("ORF")) %>%
      select(ORF,Total.TPM,pSup,SP,est_counts.Total,est_counts.Sup,est_counts.Pellet,
             Mixing_ratio.Sup,counts_adj.Sup,Mixing_ratio.Pellet,counts_adj.Pellet)

    # calcualte Zsup

    odds <- function(p) {
      p/(1-p)
    }
    logodds <- function(x) {
      # log odds, a shortcut
      res <- log(odds(x))
      res[!is.finite(res)] <- NA
      res
    }
    rollup <- function(x, flds, name, ...) {
      y <- zoo(x[,flds])
      res <- as.data.frame(rollapply(zoo(x[,flds]),...))
      colnames(res) <- paste(name,flds,sep='.')
      as_tibble(res) %>% mutate(ORF=x$ORF)
    }
    binwidth <- 100
    min.binwidth <- 5 # For edges

    # Sort by length so bins follow length
    df_mean.lo.psup <- df_pSup %>%
      left_join(gene_labels,by="ORF") %>%
      filter(classification=="Verified") %>%
      arrange(LengthTxEst) %>%
      mutate(lo.psup = logodds(pSup)) %>%
      ungroup %>%
      nest %>%
      mutate(data=map(data,function(df){
        roll.mean <- rollup(df, c('lo.psup'), "mean",width=binwidth, FUN=mean, by.column=TRUE, na.rm=TRUE, fill=NA, partial=min.binwidth)
        roll.sd <- rollup(df, c('lo.psup'), "sd",width=binwidth, FUN=sd, by.column=TRUE, na.rm=TRUE, fill=NA, partial=min.binwidth)
        return(roll.mean %>% full_join(roll.sd,by="ORF"))
      })) %>% unnest(data)
    df_Zsup <- df_pSup %>%
      left_join(df_mean.lo.psup,by=c("ORF")) %>%
      mutate(lo.psup=logodds(pSup),
             Zsup = (lo.psup-mean.lo.psup)/sd.lo.psup) %>%
      select(ORF,Total.TPM,pSup,SP,Zsup,lo.psup,mean.lo.psup,sd.lo.psup,
             est_counts.Total,est_counts.Sup,est_counts.Pellet,
             Mixing_ratio.Sup,counts_adj.Sup,Mixing_ratio.Pellet,counts_adj.Pellet)

    return(df_Zsup)
  })) %>%
  unnest(data) %>%
  left_join(df_samples_byLysate,by="Lysate_sample")

write_tsv(df_Zsup,file.path(github_dir,"data_processed","sedseq_out.tsv.gz"))

df_pSup_pombe <- df_samples %>%
  filter(Treatment_group=="spombe_spikeIn") %>%
  left_join(ratios_join,by=c("Dataset","Lysate_sample","Fraction")) %>%
  select(Dataset,Lysate_sample,Fraction,Kallisto_file,Mixing_ratio) %>%
  group_by(Dataset,Lysate_sample) %>%
  nest %>%
  mutate(data=map(data,function(df){
    df_counts <- df$Kallisto_file %>% map(function(file) {read_kallisto(file)}) %>%
      bind_rows %>% left_join(df_samples,by="Kallisto_file") %>%
      left_join(df,by=c("Kallisto_file","Fraction"))
    df_tot <- df_counts %>% filter(Fraction=="Total") %>%
      mutate(Total.TPM = TPM,
             est_counts.Total=est_counts) %>%
      select(ORF,Total.TPM,est_counts.Total)
    df_pSup <- df_counts %>% filter(Fraction!="Total") %>%
      select(ORF,Fraction,est_counts,Mixing_ratio) %>%
      mutate(counts_adj = est_counts*Mixing_ratio) %>%
      pivot_wider(names_from=Fraction,values_from=c(est_counts,counts_adj,Mixing_ratio),names_sep=".") %>%
      mutate(mixSum=counts_adj.Sup+counts_adj.Pellet,
             pSup=counts_adj.Sup/mixSum,
             SP = counts_adj.Sup/counts_adj.Pellet) %>%
      left_join(df_tot,
                by=c("ORF")) %>%
      select(ORF,Total.TPM,pSup,SP,est_counts.Total,est_counts.Sup,est_counts.Pellet,
             Mixing_ratio.Sup,counts_adj.Sup,Mixing_ratio.Pellet,counts_adj.Pellet)
  })) %>% unnest(data) %>%
  write_tsv(file.path(github_dir,"data_processed","pSup_spombe_spikeIn.tsv.gz"))

df_Zsup_filt <- df_Zsup %>%
  left_join(gene_labels,by="ORF") %>%
  filter(classification=="Verified") %>%
  group_by(Lysate_sample) %>%
  mutate(Total.TPM = Total.TPM/sum(Total.TPM,na.rm=T)*1e6) %>%
  ungroup %>%
  filter(Total.TPM>1) %>%
  group_by(ORF,LengthTxEst,Strain,Treatment_group,Treatment,Treated,Temperature) %>%
  filter(length(ORF)==2|(Treatment_group=="hairpinReporters")|
           Treatment_group=="CHX"&Temperature=="42C") %>%
  summarise(Total.TPM.mean = exp(mean(log(Total.TPM))),
            pSup.mean = mean(pSup),
            SP.mean = exp(mean(log(SP))),
            Zsup.mean = mean(Zsup)) %>%
  write_tsv(file.path(github_dir,"data_processed","sedseq_filt_mean.tsv.gz"))


