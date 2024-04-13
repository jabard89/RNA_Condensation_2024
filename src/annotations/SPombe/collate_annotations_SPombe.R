library(tidyverse)
gene_table <- read_tsv("Schizosaccharomyces_pombe_230503_gene_IDs_names_products.tsv") %>%
  rename("ORF" = "gene_systematic_id")
cds_length <- read_tsv("Schizosaccharomyces_pombe_cds_length.tsv",comment="#") %>%
  mutate(Feature="CDS")
UTR5_length <- read_tsv("Schizosaccharomyces_pombe_230503_5UTR_length.tsv",comment="#") %>%
  mutate(Feature="UTR5")
UTR3_length <- read_tsv("Schizosaccharomyces_pombe_230503_3UTR_length.tsv",comment="#") %>%
  mutate(Feature="UTR3")

length_table <- bind_rows(cds_length,UTR5_length,UTR3_length) %>%
  pivot_wider(names_from="Feature",values_from="Length",names_prefix="length.")

median_lengths <- length_table %>%
  left_join(gene_table,by="ORF") %>%
  filter(gene_type=="protein coding gene") %>%
  summarise(length.UTR5.median = median(length.UTR5,na.rm=T),
            length.UTR3.median = median(length.UTR3,na.rm=T))

# use the median UTR5 and UTR3 length as an estimate
gene_table_out <- gene_table %>%
  left_join(length_table,by="ORF") %>%
  mutate(length.UTR5.Est = if_else(!is.na(length.UTR5),length.UTR5,
                                   if_else(gene_type=="protein coding gene",median_lengths$length.UTR5.median,0)),
         length.UTR3.Est = if_else(!is.na(length.UTR3),length.UTR3,
                                   if_else(gene_type=="protein coding gene",median_lengths$length.UTR3.median,0)),
         LengthTxEst = length.CDS + length.UTR5.Est + length.UTR3.Est)

write_tsv(gene_table_out,"Schizosaccharomyces_pombe_230503_gene_IDs_names_products_lengths.tsv")
