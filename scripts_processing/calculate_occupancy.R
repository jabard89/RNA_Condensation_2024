# Normalize PolySeq by spike-in
# Jared Bard
# updated 240414
library(tidyverse)
library(conflicted)
library(here)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
github_dir <- here()
df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))) %>%
  mutate(Dataset = paste0("snake_",Sequencing_date))
df_samples_byLysate <- read_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))

Scer_labels <- read_tsv(file.path(github_dir,"src/annotations/labeled_genes_scer.tsv")) %>%
  mutate(Species="Scer")
Spombe_labels <- read_tsv(file.path(github_dir,"src/annotations/SPombe/Schizosaccharomyces_pombe_230503_gene_IDs_names_products_lengths.tsv")) %>%
  filter(gene_type=="protein coding gene") %>%
  mutate(classification="Verified") %>%
  mutate(Species="Spombe")
gene_labels_comb <- Scer_labels %>%
  bind_rows(Spombe_labels %>% select(ORF,Species,classification))
df_PolySeq_counts <- read_tsv(file.path(github_dir,"data_raw","all_kallisto_counts.tsv.gz")) %>%
  left_join(df_samples,by=c("Kallisto_file","RNA_sample")) %>%
  dplyr::filter(Experiment=="PolySeq")
df_PolySeq_TPM <- df_PolySeq_counts %>%
  filter(Fraction=="Total") %>%
  mutate("Total.TPM"=TPM) %>%
  select(ORF,Dataset,Lysate_sample,Total.TPM)

df_PolySeq_Pombe <- df_PolySeq_counts %>%
  filter(Dataset!="snake_210512_210520") %>%
  left_join(gene_labels_comb %>% select(ORF,Species),by="ORF") %>%
  filter(Species=="Spombe") %>%
  select(ORF,Dataset,RNA_sample,est_counts) %>%
  filter(est_counts>100) %>%
  group_by(Dataset) %>%
  nest %>%
  mutate(data=map(data,function(df){
    df %>% pivot_wider(names_from=RNA_sample,values_from=est_counts) %>%
      drop_na %>%
      pivot_longer(cols=c(-ORF),names_to="RNA_sample",values_to="est_counts") %>%
      group_by(ORF) %>%
      mutate(norm.counts = est_counts/median(est_counts)) %>%
      group_by(RNA_sample) %>%
      summarise(Correction.factor.med = median(norm.counts),
                Correction.factor.sd = sd(norm.counts),
                percent_error=Correction.factor.sd/Correction.factor.med)
  })) %>%
  unnest(data)
write_tsv(df_PolySeq_Pombe %>%
            left_join(df_samples,by=c("Dataset","RNA_sample")),
          file.path(github_dir,"data_processed","mixingratios_Pombespike_Polyseq.tsv"))
df_PolySeq_normPombe <- df_PolySeq_counts %>%
  left_join(Scer_labels,by="ORF") %>%
  filter(Species=="Scer") %>%
  filter(Dataset!="snake_210512_210520") %>%
  left_join(df_PolySeq_Pombe,by=c("Dataset","RNA_sample")) %>%
  mutate(Counts.norm=est_counts/Correction.factor.med) %>%
  filter(est_counts>10) %>%
  select(ORF,Dataset,Strain,Lysate_sample,Treatment_group,Treatment,Treated,Temperature,Time,Biorep,Fraction,Counts.norm,est_counts) %>%
  pivot_wider(names_from=Fraction,values_from=c(Counts.norm,est_counts),names_sep="_") %>%
  mutate(Counts.norm_Mono=if_else(is.na(Counts.norm_Mono),0,Counts.norm_Mono),
         est_counts_Mono=if_else(is.na(est_counts_Mono),0,est_counts_Mono)) %>%
  drop_na

df_PolySeq_normSpike <- df_PolySeq_counts %>%
  filter(Dataset=="snake_210512_210520") %>%
  group_by(RNA_sample) %>%
  mutate(Counts.norm=est_counts/est_counts[ORF=="SSA2_Nanoluc_TPI1"]) %>%
  ungroup %>%
  select(ORF,Dataset,Strain,Lysate_sample,Treatment_group,Treatment,Treated,Temperature,Time,Biorep,Fraction,Counts.norm,est_counts) %>%
  pivot_wider(names_from=Fraction,values_from=c(Counts.norm,est_counts))
rm(df_PolySeq_counts)

df_PolySeq_Occ <- bind_rows(df_PolySeq_normPombe,df_PolySeq_normSpike)
rm(df_PolySeq_normPombe,df_PolySeq_normSpike)
df_PolySeq_Occ <- df_PolySeq_Occ %>%
  left_join(Scer_labels,by="ORF") %>%
  left_join(df_PolySeq_TPM,by=c("ORF","Dataset","Lysate_sample")) %>%
  mutate(Counts.norm_Mono=if_else(is.na(Counts.norm_Mono),0,Counts.norm_Mono)) %>%
  mutate(Occ = (Counts.norm_Mono+Counts.norm_Poly)/(Counts.norm_Free + Counts.norm_Mono + Counts.norm_Poly),
         Occ_odds = (Counts.norm_Mono+Counts.norm_Poly)/(Counts.norm_Free)) %>%
  select(ORF,Dataset,Lysate_sample,
         Counts.norm_Total,Counts.norm_Free,Counts.norm_Mono,Counts.norm_Poly,
         est_counts_Total,est_counts_Free,est_counts_Mono,est_counts_Poly,
         Occ,Occ_odds,Total.TPM) %>%
  left_join(df_samples_byLysate,by="Lysate_sample")
write_tsv(df_PolySeq_Occ,file.path(github_dir,"data_processed","PolysomeSeq_all_normSpike.tsv.gz"))

# ggplot(df_PolySeq_Occ %>%
#          mutate(total_counts = est_counts_Free+est_counts_Mono+est_counts_Poly),
#        aes(x=total_counts,y=Total.TPM,
#            color=paste0(Temperature,".",Treatment)))+
#   geom_point()+
#   scale_x_log10nice()+scale_y_log10nice()+
#   facet_grid(Treatment_group~Biorep)+
#   geom_hline(yintercept=1)+
#   geom_vline(xintercept=50)+
#   coord_cartesian(xlim=c(10,1e5),ylim=c(0.1,1e5))

df_PolySeq_filt_mean <- df_PolySeq_Occ %>%
  left_join(Scer_labels,by="ORF") %>% filter(classification=="Verified") %>%
  group_by(Lysate_sample) %>%
  mutate(Total.TPM = Total.TPM/sum(Total.TPM,na.rm=T)*1e6) %>%
  ungroup %>%
  mutate(total_counts = est_counts_Free+est_counts_Mono+est_counts_Poly) %>%
  group_by(ORF,Strain,Treatment_group,Treatment,Treated,Temperature,Time) %>%
  filter(total_counts>10) %>%
  summarise(Total.TPM.Poly.mean = exp(mean(log(Total.TPM))),
            Occ.mean = mean(Occ),
            Occ.odds.mean = exp(mean(log(Occ_odds)))) %>%
  write_tsv(file.path(github_dir,"data_processed","PolySeq_minfilt_mean.tsv.gz"))
