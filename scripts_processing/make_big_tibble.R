library(tidyverse)
library(here)

github_dir <- here()
source(file.path(github_dir,"scripts_plotting/main_stylesheet_230605.R"))
gene_labels <- gene_labels <- read_tsv(file.path(github_dir,"src/annotations/labeled_genes_scer.tsv")) %>%
  dplyr::mutate(label=case_when(label=="RiBi"~"ribosome biogenesis",
                                label=="translation factors"~"other",
                                label=="RP"~"other",
                                label=="glycolysis"~"other",
                                label=="HSF1"~"HSF1 targets",
                                label=="MSN2/4"~"MSN2/4 targets",
                                TRUE~label)) %>%
  dplyr::mutate(label=factor(label,levels=label.levels)) %>%
  dplyr::mutate(length.trans = LengthTxEst) %>%
  select(ORF,gene,classification,label,LengthTxEst,length.protein)

df_samples <- read_csv(file.path(github_dir,"scripts_processing/fastq_processing/fastq_table.csv")) %>%
  mutate(Kallisto_file = file.path(github_dir,"data_raw/kallisto_quant/tsv",paste0(RNA_sample,"~abundance.tsv"))) %>%
  mutate(Dataset = paste0("snake_",Sequencing_date))
df_samples_byLysate <- read_tsv(file.path(github_dir,"data_raw","RNAseq_samples_byLysate.tsv"))
df_Occ_mean <- read_tsv(file.path(github_dir,"data_processed","PolySeq_minfilt_mean.tsv.gz"))

df_Zsup_mean <- read_tsv(file.path(github_dir,"data_processed","sedseq_filt_mean.tsv.gz")) %>%
  left_join(df_samples_byLysate %>% select(Treatment_group,Treatment,Temperature,Time) %>% unique,by=c("Treatment","Treatment_group","Temperature"))
df_stress_samples <- read_csv(file.path(github_dir,"data_raw/stress_labels.csv")) %>%
  mutate(Stress=factor(Stress,levels=c("none","Heat Shock","Azide","Ethanol","DTT")),
         Stress_label=factor(Stress_label,levels=c("none","30C","42C","46C",
                                                   "mock","Azide 0.5%","Azide 0.8%",
                                                   "Ethanol 5%","Ethanol 7.5%","Ethanol 10%","Ethanol 15%",
                                                   "DTT 10mM")))
df_counts <- read_tsv((file.path(github_dir,"data_raw","all_kallisto_counts.tsv.gz")),
                      comment="#") %>%
  left_join(df_samples,by=c("Kallisto_file","RNA_sample")) %>%
  left_join(df_stress_samples,by=c("Temperature","Treatment_group","Treatment"))

df_Occ_mean_minfilt <- read_tsv(file.path(github_dir,"data_processed","PolySeq_minfilt_mean.tsv.gz")) %>%
  left_join(df_stress_samples,by=c("Temperature","Treatment_group","Treatment")) %>%
  filter(!is.na(Control)) %>%
  group_by(Treatment_group,ORF) %>%
  filter(length(ORF[Control==TRUE])==1) %>%
  mutate(Occ.odds.FC = Occ.odds.mean/Occ.odds.mean[Control==TRUE])
df_poly_mean <- read_tsv(file.path(github_dir,"data_processed","deseq_polysomes_mean.tsv.gz")) %>%
  left_join(df_stress_samples,by=c("Temperature","Treatment_group","Treatment")) %>%
  filter(!is.na(Stress_group)) %>%
  left_join(df_samples_byLysate %>% ungroup %>% select(Treatment_group,Time) %>% unique,by="Treatment_group")

polyseq_azide_deseq <- read_tsv(file.path(github_dir,"data_processed/deseq_results","230508_PolySeq_0.8%Azide_deseq2.tsv")) %>%
  rename("TPM.FC.Poly"="FC.vs.mock") %>%
  mutate(Temperature="30C")
polyseq_etoh_deseq <- read_tsv(file.path(github_dir,"data_processed/deseq_results","230508_PolySeq_7.5%EtOH_deseq2.tsv")) %>%
  rename("TPM.FC.Poly"="FC.vs.mock") %>%
  mutate(Temperature="30C")
poly10_deseq <- read_tsv(file.path(github_dir,"data_processed/deseq_results","230508_PolySeq_10min_deseq2.tsv")) %>%
  rename("TPM.FC.Poly"="FC.vs.30C") %>%
  mutate(Treatment="none")
sedseq_azide3_deseq <- read_tsv(file.path(github_dir,"data_processed/deseq_results","azide3_deseq2_230510.tsv")) %>%
  rename("TPM.FC.SedSeq"="FC.vs.mock") %>%
  mutate(Temperature="30C")
sedseq_milderEtOH2_deseq <- read_tsv(file.path(github_dir,"data_processed/deseq_results","milderEtOH2_deseq2_230510.tsv")) %>%
  rename("TPM.FC.SedSeq"="FC.vs.mock") %>%
  mutate(Temperature="30C")
sedseq_compTemp_bygroup <- read_tsv(file.path(github_dir,"data_processed/deseq_results","TSP_DESeq_compTemp_byGroup_deseq2_230510.tsv")) %>%
  rename("TPM.FC.SedSeq"="FC.vs.30C")
sedseq_azide0dot5_deseq <-  read_tsv(file.path(github_dir,"data_processed/deseq_results","azide0dot5_deseq2_230525.tsv")) %>%
  rename("TPM.FC.SedSeq"="FC.vs.mock") %>%
  mutate(Temperature="30C")
deseq_polyseq <- bind_rows(polyseq_azide_deseq,polyseq_etoh_deseq,poly10_deseq) %>%
  select(-padj)
deseq_sedseq <- bind_rows(sedseq_azide3_deseq,sedseq_milderEtOH2_deseq,sedseq_compTemp_bygroup,sedseq_azide0dot5_deseq) %>%
  select(-padj)


df_Zsup_Poly_minfilt <- df_Zsup_mean  %>%
  full_join(df_stress_samples,by=c("Temperature","Treatment_group","Treatment")) %>%
  mutate(Control=if_else(Stress_group=="hairpinReporters"&Strain=="yHG010",TRUE,Control)) %>%
  filter(!is.na(Stress_group)) %>%
  group_by(ORF,Strain,Treatment_group,Stress,Stress_group,Stress_label,Control,Temperature,Treatment) %>%
  summarise(TPM.SedSeq = exp(mean(log(Total.TPM.mean))),
            SP.mean  = exp(mean(log(SP.mean))),
            pSup.mean = SP.mean/(1+SP.mean),
            Zsup.mean = mean(Zsup.mean),
            esc.mean = mean(esc.mean),
            sed.mean = mean(sed.mean)) %>%
  group_by(ORF,Stress,Stress_group) %>%
  filter(length(ORF[Control==TRUE])==1) %>%
  mutate(SP.FC = SP.mean/SP.mean[Control==TRUE],
         deltaZ = Zsup.mean-Zsup.mean[Control==TRUE]) %>%
  left_join(deseq_sedseq,by=c("ORF","Treatment_group","Temperature","Treatment")) %>%
  ungroup %>%
  select(-Treatment_group) %>%
  full_join(df_poly_mean %>% ungroup %>% filter(Treatment_group%in%c("PolySeq_20min","PolySeq_10min",
                                                                     "PolySeq_7.5%EtOH","PolySeq_0.8%Azide")) %>%
              rename("TPM.Total.Poly"="TPM.Total") %>%
              select(ORF,Stress,Stress_group,Stress_label,Treatment,Temperature,Time,Control,
                     TPM.Total.Poly,Rel_Occ,FC.vs.ctrl_RNA,FC.vs.ctrl_Rib,FC.vs.ctrl_TE,sig_RNA,sig_Rib,sig_TE),
            by=c("ORF","Stress","Stress_group","Stress_label","Control","Temperature","Treatment")) %>%
  full_join(df_Occ_mean_minfilt %>% ungroup %>% filter(Treatment_group%in%c("PolySeq_20min","PolySeq_10min",
                                                                            "PolySeq_7.5%EtOH","PolySeq_0.8%Azide")) %>%
              select(ORF,Stress,Stress_group,Stress_label,Treatment,Temperature,Time,Control,
                     Occ.mean,Occ.odds.mean,Occ.odds.FC),
            by=c("ORF","Stress","Stress_group","Stress_label","Control","Temperature","Treatment","Time")) %>%
  left_join(gene_labels,by=c("ORF")) %>%
  mutate(Stress_group=case_when(
    Stress_group=="Heat Shock"&Time=="20min"~"HS.20minPoly",
    TRUE~Stress_group
  )) %>%
  arrange(label) %>%
  filter(SP.mean>0,
         SP.FC>0,
         !is.infinite(SP.FC))
write_tsv(df_Zsup_Poly_minfilt,file.path(github_dir,"data_processed","df_Zsup_Poly_minfilt.tsv.gz"))
