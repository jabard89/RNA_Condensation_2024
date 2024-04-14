library(tidyverse)

src_dir <- "/Users/Caitlin/Dropbox (Drummond Lab)/smFISH_data/fishquant"
dates <- c("20230810", "20230816", "20230823", "20230825")
transcripts <- c("ADD66", "SSA4", "SSB1", "HSP104")
conditions <- c("30C", "42C", "46C")

sample_table_list <- list()

for (date in dates) {
  date_dir <- file.path(src_dir, date)
  
  for (trans in list.dirs(date_dir, recursive = FALSE)) {
    trans_basename <- basename(trans)
    if (trans_basename %in% transcripts) {
      trans_dir <- file.path(date_dir, trans_basename)
      
      for (cond in list.dirs(trans_dir, recursive = FALSE)) {
        cond_basename <- basename(cond)
        if (cond_basename %in% conditions) {
          cond_dir <- file.path(trans_dir, cond_basename)
          out_row <- tibble("Date" = date, "Transcript" = trans_basename, "Condition" = cond_basename, "Condition.dir" = cond_dir)
          sample_table_list <- c(sample_table_list, list(out_row))
        }
      }
    }
  }
}

sample_table <- bind_rows(sample_table_list)

df_rna <- pmap(sample_table,function(Date,Transcript,Condition,Condition.dir) {
  file_name = paste0(Date,"~",Transcript,"~",Condition,"_rna_df.tsv.gz")
  file <- file.path(Condition.dir,file_name)
  read_tsv(file)  %>%
    mutate(Image_id = as.character(Image_id))
}) %>% bind_rows() %>%
  pivot_longer(cols=starts_with("RNA_mean"),
               names_pattern = "(.*)_mean_(.*)",
               names_to=c("Key","Channel"),
               values_to="ROI_mean") %>%
  rename("loc_x"="RNA_loc_x",
         "loc_y"="RNA_loc_y")
write_tsv(df_rna,file.path(src_dir,"230930_analyze_smfish_RNA.tsv.gz"))
df_random <- pmap(sample_table,function(Date,Transcript,Condition,Condition.dir) {
  file_name = paste0(Date,"~",Transcript,"~",Condition,"_random_df.tsv.gz")
  file <- file.path(Condition.dir,file_name)
  read_tsv(file)  %>%
    mutate(Image_id = as.character(Image_id))
}) %>% bind_rows %>%
  pivot_longer(cols=starts_with("random_mean"),
               names_pattern = "(.*)_mean_(.*)",
               names_to=c("Key","Channel"),
               values_to="ROI_mean") %>%
  rename("loc_x"="random_loc_x",
         "loc_y"="random_loc_y")
write_tsv(df_random,file.path(src_dir,"230930_analyze_smFISH_random.tsv.gz"))