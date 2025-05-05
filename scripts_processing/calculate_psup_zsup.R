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

github_dir <- here::here()
source(file.path(github_dir,"scripts_processing/utilityFunctions.R"))
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

    # calculate Zsup

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

# Plan:
# 1. Manually identify "control" lysate
# 2. Calculate rolling mean using windowfunction
# 3. calculate sed and esc per rep
df_Zsup_filt <- df_Zsup %>% ungroup %>%
  left_join(gene_labels %>% select(ORF,gene,LengthTxEst,classification),by="ORF") %>%
  filter(classification=="Verified") %>%
  group_by(Lysate_sample) %>%
  mutate(Total.TPM = Total.TPM/sum(Total.TPM,na.rm=T)*1e6) %>%
  filter(Total.TPM > 1,
         pSup > 0,
         !is.infinite(pSup)) %>% ungroup

lengths = sort(unique(as.integer((df_Zsup_filt %>% pull(LengthTxEst)))))
df_pSup_windowmean <- df_Zsup_filt %>%
  group_by(Lysate_sample) %>%
  nest %>%
  mutate(data = map(data, function(df) {
    windowfunction(df,"lo.psup", "LengthTxEst", "Treatment",
                   n=lengths, minwidth_fraction=0.05, logx=TRUE, minn=1, na.rm=T) %>%
      mutate("pSup.lysate.windowmean" = invlogodds(lo.psup)) %>% select(LengthTxEst,pSup.lysate.windowmean)
    })) %>%
  unnest(c(data))

df_ctrl <- df_Zsup_filt %>% ungroup %>%
  left_join(df_pSup_windowmean, by = c("Lysate_sample", "LengthTxEst")) %>%
  group_by(Treatment_group,ORF) %>%
  mutate(pSup.ctrl.mean = invlogodds(mean(logodds(pSup[Control == TRUE]))),
         pSup.ctrl.windowmean.mean = invlogodds(mean(logodds(pSup.lysate.windowmean[Control == TRUE])))) %>%
  group_by(Treatment_group) %>%
  mutate(pSup.ctrl.sd = sd(logodds(pSup.ctrl.mean) - logodds(pSup.ctrl.windowmean.mean), na.rm=T)) %>%
  select(ORF,Treatment_group,pSup.ctrl.mean,pSup.ctrl.windowmean.mean,pSup.ctrl.sd) %>%
  unique()

df_Zsup_sed_esc <- df_Zsup %>% ungroup %>%
  left_join(gene_labels %>% select(ORF,gene,LengthTxEst,classification),by="ORF") %>%
  left_join(df_pSup_windowmean, by = c("Lysate_sample", "LengthTxEst")) %>%
  left_join(df_ctrl, by = c("ORF","Treatment_group")) %>%
  group_by(Lysate_sample) %>%
  mutate(pSup.lysate.sd = sd(logodds(pSup) - logodds(pSup.lysate.windowmean), na.rm=T)) %>%
  ungroup %>%
  mutate(sed = (logodds(pSup.ctrl.mean) - logodds(pSup)) / pSup.ctrl.sd,
         esc = ((logodds(pSup) - logodds(pSup.lysate.windowmean)) - (logodds(pSup.ctrl.mean) - logodds(pSup.ctrl.windowmean.mean))) / pSup.ctrl.sd,
         RelSed = (logodds(pSup) - logodds(pSup.lysate.windowmean)) / pSup.lysate.sd,
         pSup.ctrl.window.mean=pSup.ctrl.windowmean.mean,
         pSup.treatment.window.mean=pSup.lysate.windowmean,
         lopSup.treatment.sd=pSup.lysate.sd,
         lopSup.ctrl.sd=pSup.ctrl.sd)


write_tsv(df_Zsup_sed_esc,file.path(github_dir,"data_processed","sedseq_out.tsv.gz"))

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

df_Zsup_sed_esc_filt_mean <- df_Zsup_sed_esc %>%
  filter(classification=="Verified") %>%
  group_by(Lysate_sample) %>%
  mutate(Total.TPM = Total.TPM/sum(Total.TPM,na.rm=T)*1e6) %>%
  ungroup %>%
  filter(Total.TPM>1, pSup>0, !is.infinite(pSup)) %>%
  group_by(ORF,LengthTxEst,Strain,Treatment_group,Treatment,Treated,Temperature) %>%
  filter(length(ORF)==2|(Treatment_group=="hairpinReporters")|
           Treatment_group=="CHX"&Temperature=="42C") %>%
  summarise(Total.TPM.mean = exp(mean(log(Total.TPM))),
            pSup.mean = invlogodds(mean(logodds(pSup))),
            SP.mean = exp(mean(log(SP))),
            Zsup.mean = mean(Zsup),
            sed.mean = mean(sed),
            esc.mean = mean(esc),
            RelSed.mean = mean(RelSed),
            pSup.ctrl.window.mean.mean=invlogodds(mean(logodds(pSup.ctrl.windowmean.mean))),
            pSup.treatment.window.mean.mean=invlogodds(mean(logodds(pSup.lysate.windowmean))),
            lopSup.treatment.sd.mean=mean(pSup.lysate.sd),
            lopSup.ctrl.sd.mean=mean(pSup.ctrl.sd)
  ) %>%
  write_tsv(file.path(github_dir,"data_processed","sedseq_filt_mean.tsv.gz"))


