# Normalize PolySeq counts by spike-in first, then combine polysome and monosome counts
# then use deseq2 to calculate fold changes in "TE"
# Jared Bard
# 230505
library(tidyverse)
library(DESeq2)
library(conflicted)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
github_dir <- here()
df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))) %>%
  mutate(Dataset = paste0("snake_",Sequencing_date)) %>%
  mutate(Temperature=factor(Temperature,levels=c("30C","37C","40C","42C","46C")),
         Treatment=factor(Treatment,levels=c("none","mock","Azide_0.8%_pH6.8","EtOH_7.5%")))
df_samples_byLysate <- read_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))
gene_labels <- read_tsv(file.path(github_dir,"src/annotations/labeled_genes_scer.tsv")) 

df_PolySeq_counts <- read_tsv(file.path(github_dir,"data_raw","all_kallisto_counts.tsv.gz")) %>%
  left_join(df_samples,by=c("Kallisto_file","RNA_sample")) %>%
  dplyr::filter(Experiment=="PolySeq")

# load correction factors to allow for mono and poly to be combined

df_Poly_spike <- read_tsv(file.path(github_dir,"data_processed","PolysomeSeq_all_normSpike.tsv.gz"))
df_pombe_mixingratios <- read_tsv(file.path(github_dir,"data_processed","mixingratios_Pombespike_Polyseq.tsv")) %>%
  filter(Dataset!="snake_221019")%>%
  select(RNA_sample,Correction.factor.med)
df_PolySeq_normSpike <- df_PolySeq_counts %>%
  filter(Dataset=="snake_210512_210520") %>%
  filter(ORF=="SSA2_Nanoluc_TPI1") %>%
  group_by(Lysate_sample) %>%
  mutate(Correction.factor.med = est_counts/est_counts[Fraction=="Poly"]) %>%
  ungroup %>%
  select(RNA_sample,Correction.factor.med)
corr_factors <- bind_rows(df_pombe_mixingratios,df_PolySeq_normSpike) %>%
  left_join(df_samples,by="RNA_sample") %>%
  filter(Dataset!="snake_221019")%>%
  group_by(Lysate_sample) %>%
  mutate(Correction.factor.med = Correction.factor.med/Correction.factor.med[Fraction=="Poly"],
         Correction.factor.med = if_else(Fraction=="Total",1,Correction.factor.med)) %>%
  ungroup %>%
  select(RNA_sample,Correction.factor.med)
df_PolySeq_corr <- df_PolySeq_counts %>%
  left_join(gene_labels,by="ORF") %>%
  filter(classification=="Verified") %>%
  left_join(corr_factors,by="RNA_sample") %>%
  select(ORF,Dataset,Lysate_sample,Fraction,est_counts,Correction.factor.med,length.protein) %>%
  filter(est_counts>0) %>%
  group_by(Dataset,Lysate_sample,Fraction) %>%
  mutate(counts.corr = est_counts/Correction.factor.med) %>%
  select(-Correction.factor.med,-est_counts) %>%
  pivot_wider(names_from=Fraction,values_from=counts.corr,names_prefix="Counts.") %>%
  replace_na(list(Counts.Mono=0)) %>%
  mutate(Counts.Rib = Counts.Mono+Counts.Poly) %>%
  pivot_longer(starts_with("Counts"),names_to="Fraction",values_to="Counts.corr",names_prefix = "Counts.") %>%
  group_by(Lysate_sample,Fraction) %>%
    mutate(RPKB = Counts.corr/(length.protein*3/1e3),
           TPM = RPKB/sum(RPKB,na.rm=T)*1e6) %>%
  select(Lysate_sample,ORF,Fraction,TPM,Counts.corr) %>%
  pivot_wider(names_from=Fraction,values_from=c("TPM","Counts.corr"))
write_tsv(df_PolySeq_corr,file.path(github_dir,"data_processed","Polysome_counts_raw.tsv.gz"))

df_PolySeq_counts_corr_long <- df_PolySeq_corr %>%
  select(ORF,Lysate_sample,starts_with("Counts")) %>%
  pivot_longer(cols=starts_with("Counts"),names_to="Fraction",values_to="Counts.corr",names_prefix = "Counts.corr_")

poly_samples <- df_PolySeq_counts_corr_long %>%
  select(Lysate_sample,Fraction) %>% unique %>%
  left_join(df_samples_byLysate,by="Lysate_sample")

calculate_deseq <- function(t_group,contrast_var,baseline){
  samples <- bind_rows(poly_samples %>% mutate(Fraction="Total"),
                      poly_samples %>% mutate(Fraction="Rib"))
  contrast_var_sym <- sym(contrast_var)
  all_con_vars <- poly_samples %>%
    dplyr::filter(Treatment_group==!!t_group) %>%
    dplyr::filter(Fraction%in%c("Total","Rib")) %>%
    unite(col="Sample",Lysate_sample,Fraction) %>%
    mutate("Contrast_variable"=factor({{contrast_var_sym}})) %>%
    pull(Contrast_variable) %>% unique %>% as.character
  non_baseline_values <- all_con_vars[all_con_vars!=baseline]
  
  cts_all <- df_PolySeq_counts_corr_long %>%
    select(ORF,Lysate_sample,Fraction,Counts.corr) %>%
    unite(col="Sample",Lysate_sample,Fraction) %>%
    mutate(Counts.corr = round(Counts.corr)) %>%
    filter(!is.na(Counts.corr),Counts.corr>0) %>%
    pivot_wider(names_from=Sample,values_from=Counts.corr) %>%
    column_to_rownames(var="ORF") %>% as.matrix(.)
  
 cts_all[is.na(cts_all)] <- 0
  
  coldata_TE <- poly_samples %>%
    dplyr::filter(Treatment_group==!!t_group) %>%
    dplyr::filter(Fraction%in%c("Total","Rib")) %>%
    unite(col="Sample",Lysate_sample,Fraction,remove=F) %>%
    column_to_rownames(var="Sample") %>%
    mutate(Fraction=factor(Fraction,levels=c("Total","Rib"))) %>%
    mutate("Contrast_variable"=factor({{contrast_var_sym}}))
  
  cts_TE <- cts_all[,rownames(coldata_TE)]
  dds_TE <- DESeqDataSetFromMatrix(countData=cts_TE,
                                colData=coldata_TE,
                                design = ~ Fraction + Contrast_variable + Fraction:Contrast_variable) %>%
    DESeq(.) %>% estimateSizeFactors(.)
  df_results_TE <- non_baseline_values %>% map(function(numerator){
    FC_col_name <- paste0("FC.vs.",{{baseline}})
    df <- results(dds_TE,tidy=T,alpha=0.05,name=paste0("FractionRib.Contrast_variable",numerator))
    mutate(df,{{contrast_var}} := numerator) %>%
      dplyr::rename("ORF"="row") %>%
      mutate(FC.vs.ctrl := 2^log2FoldChange) %>%
      select(ORF,FC.vs.ctrl,padj,{{contrast_var}}) %>%
      mutate(type="TE")
  }) %>% bind_rows
  
  coldata_RNA <- poly_samples %>%
    dplyr::filter(Treatment_group==!!t_group) %>%
    dplyr::filter(Fraction%in%c("Total")) %>%
    unite(col="Sample",Lysate_sample,Fraction,remove=F) %>%
    column_to_rownames(var="Sample") %>%
    mutate("Contrast_variable"=factor({{contrast_var_sym}}))
  
  cts_RNA <- cts_all[,rownames(coldata_RNA)]
  dds_RNA <- DESeqDataSetFromMatrix(countData=cts_RNA,
                                colData=coldata_RNA,
                                design = ~ Contrast_variable) %>%
    DESeq(.) %>% estimateSizeFactors(.)
  df_results_RNA <- non_baseline_values %>% map(function(numerator){
    FC_col_name <- paste0("FC.vs.",{{baseline}})
    df <- results(dds_RNA,tidy=T,alpha=0.05,contrast=c("Contrast_variable",numerator,baseline))
    mutate(df,{{contrast_var}} := numerator) %>%
      dplyr::rename("ORF"="row") %>%
      mutate(FC.vs.ctrl := 2^log2FoldChange) %>%
      select(ORF,FC.vs.ctrl,padj,{{contrast_var}}) %>%
      mutate(type="RNA")
  }) %>% bind_rows
  
  coldata_poly <- poly_samples %>%
    dplyr::filter(Treatment_group==!!t_group) %>%
    dplyr::filter(Fraction%in%c("Rib")) %>%
    unite(col="Sample",Lysate_sample,Fraction,remove=F) %>%
    column_to_rownames(var="Sample") %>%
    mutate("Contrast_variable"=factor({{contrast_var_sym}}))
  
  cts_poly <- cts_all[,rownames(coldata_poly)]
  dds_poly <- DESeqDataSetFromMatrix(countData=cts_poly,
                                    colData=coldata_poly,
                                    design = ~ Contrast_variable) %>%
    DESeq(.) %>% estimateSizeFactors(.)
  df_results_poly <- non_baseline_values %>% map(function(numerator){
    FC_col_name <- paste0("FC.vs.",{{baseline}})
    df <- results(dds_poly,tidy=T,alpha=0.05,contrast=c("Contrast_variable",numerator,baseline))
    mutate(df,{{contrast_var}} := numerator) %>%
      dplyr::rename("ORF"="row") %>%
      mutate(FC.vs.ctrl := 2^log2FoldChange) %>%
      select(ORF,FC.vs.ctrl,padj,{{contrast_var}}) %>%
      mutate(type="Rib")
  }) %>% bind_rows
  
  bind_rows(df_results_RNA,df_results_poly,df_results_TE)
}

deseq_10min_raw <- calculate_deseq("PolySeq_10min","Temperature","30C") %>%
  mutate(Treatment_group="PolySeq_10min") 
deseq_10min <- deseq_10min_raw %>%
      mutate(sig = if_else(padj<0.05,TRUE,FALSE)) %>%
  select(-padj) %>%
  pivot_wider(names_from=type,values_from=c("FC.vs.ctrl","sig")) %>%
  mutate(Treated=FALSE)

deseq_azide_raw <- calculate_deseq("PolySeq_0.8%Azide","Treated","FALSE") %>%
  mutate(Treatment_group="PolySeq_0.8%Azide")
deseq_azide <- deseq_azide_raw %>%
  mutate(sig = if_else(padj<0.05,TRUE,FALSE)) %>%
  select(-padj) %>%
  pivot_wider(names_from=type,values_from=c("FC.vs.ctrl","sig")) %>%
  mutate(Temperature="30C",Treated=TRUE)

deseq_etoh_raw <- calculate_deseq("PolySeq_7.5%EtOH","Treated","FALSE") %>%
  mutate(Treatment_group="PolySeq_7.5%EtOH")
deseq_etoh <- deseq_etoh_raw %>%
  mutate(sig = if_else(padj<0.05,TRUE,FALSE)) %>%
  select(-padj) %>%
  pivot_wider(names_from=type,values_from=c("FC.vs.ctrl","sig")) %>%
  mutate(Temperature="30C",Treated=TRUE)
deseq_all <- bind_rows(deseq_10min,deseq_azide,deseq_etoh)
df_poly_out <- df_PolySeq_corr %>%
  select(ORF,Lysate_sample,starts_with("TPM")) %>%
  pivot_longer(cols=starts_with("TPM"),names_to="Fraction",values_to="TPM",names_prefix = "TPM_") %>%
  left_join(df_samples_byLysate,by="Lysate_sample") %>%
  group_by(ORF,Treatment_group,Temperature,Treated,Treatment,Fraction) %>%
  summarise(TPM = exp(mean(log(TPM)))) %>%
  pivot_wider(names_from=Fraction,values_from=TPM,names_prefix="TPM.") %>%
  mutate(Rel_Occ = TPM.Rib/TPM.Total) %>%
  left_join(deseq_all,by=c("Treatment_group","ORF","Temperature","Treated"))
write_tsv(df_poly_out,file.path(github_dir,"data_processed","deseq_polysomes_mean.tsv.gz"))

