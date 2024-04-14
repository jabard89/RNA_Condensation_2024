# Jared Bard
# 230731
# Uses RNAFracQuant functions
library(RNAFracQuant)
library(tidyverse)
# Load necessary library

# Set path to the samplesheet file
samplesheet_file <- "RNAseq_samplesheet_with_kallisto_230510.tsv"

# Read the samplesheet
samplesheet <- read.table(samplesheet_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(Experiment=="SedSeq",Fraction%in%c("Total","Sup","Pellet"))

# Create the "Condition" column as a combination of "Dataset" and "Lysate_sample"
samplesheet$Condition <- paste(samplesheet$Dataset, samplesheet$Lysate_sample, sep = "_")

# Define a mapping from your fraction names to the expected fraction names
fraction_mapping <- c("Sup" = "Sup", "Pellet" = "Pellet", "Total" = "Tot")

# Rename the fractions in the samplesheet
samplesheet$Fraction <- fraction_mapping[samplesheet$Fraction]

# Create the "File" column as the path to the count file for each sample
samplesheet$File <- paste(samplesheet$Dataset, samplesheet$Lysate_sample, samplesheet$Fraction, "kallisto-counts-rounded.tsv", sep = "_")

# Write the prepared samplesheet to a file
write.table(samplesheet %>% filter(Condition=="snake_170622_170627_L7"|Condition=="snake_170622_170627_L11") %>%
              select(Temperature,Biorep,File,Fraction) ,
            "./kallisto_counts_collated/samplesheet_collated_RNAFracQuant_230731.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Get the paths to the abundance.tsv files from the samplesheet
abundance_files <- samplesheet %>% filter(Condition=="snake_170622_170627_L7"|Condition=="snake_170622_170627_L11")  %>% pull(Kallisto_file)

# Function to process a single abundance.tsv file
process_abundance_file <- function(file) {
  # Read the file
  abundance_data <- read.table(file, header = TRUE)
  
  # Extract the target_id and est_counts columns
  count_data <- abundance_data[, c("target_id", "est_counts")]
  count_data$est_counts <- round(count_data$est_counts)
  
  # Get the Dataset, Lysate_sample, and Fraction information for the file
  sample_info <- samplesheet[samplesheet$Kallisto_file == file, c("Dataset", "Lysate_sample", "Fraction")]
  
  # Convert the sample_info dataframe into a single string
  sample_info_string <- paste(unlist(sample_info), collapse = "_")
  # Create the new file name
  new_file_name <- paste(sample_info_string, "kallisto-counts-rounded.tsv", sep = "_")
  # Write the count data to a new file in the specified directory, without headers
  write.table(count_data, file.path("./kallisto_counts_collated", new_file_name), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# Process all abundance.tsv files
lapply(abundance_files, process_abundance_file)

wide_data <- RNAFracQuant::get_wide_Fraction(dir_in="./kallisto_counts_collated",file="./samplesheet_collated_RNAFracQuant_230731.txt") %>%
  rename("Condition"="Temperature")
filter_low_counts <- function(df, min_counts = 20) {
  # Select the 'Tot', 'Pellet', and 'Sup' columns
  selected_columns <- df[, c('Pellet', 'Sup')]
  
  # Calculate the row sums
  row_sums <- rowSums(selected_columns, na.rm = TRUE)
  
  # Filter out rows where the sum is less than min_sum
  # Also filter out ORFs that don't appear in all bioreps
  df_filtered <- df[row_sums >= min_counts, ] %>% filter(Tot>20) %>%
    group_by(Condition) %>%
    mutate(NREPS = length(unique(Biorep))) %>%
    group_by(Condition,ORF) %>%
    filter(length(ORF)==NREPS) %>%
    select(-NREPS)
  
  return(df_filtered)
}

matrix_list <- wide_data %>%
  filter_low_counts(.) %>%
  pivot_longer(cols=c(Tot,Sup,Pellet),names_to="Fraction",values_to="obs_count") %>%
  mutate(obs_count = round(obs_count)) %>%
  arrange(Condition,Fraction,ORF,Biorep) %>%
  pivot_wider(names_from=Biorep,values_from=obs_count) %>%
  group_by(Condition) %>% nest %>%
  mutate(data=map(data,~.x %>% group_by(Fraction) %>% split(.$Fraction) %>%
                    map(~.x%>%column_to_rownames(var="ORF")%>%select(-Fraction) %>%as.matrix(.))))



  # arrange(Temperature,Fraction,ORF) %>%
  # split(.$Temperature) %>%
  # map(~.x %>% split(.$Fraction))
  # group_by(Temperature) %>%
  # nest %>%
  # mutate(data = map(data,function(df){
  #   df %>% split(.$Fraction)}))
  #   %>%
  #     map(~ .x %>%
  #           split(.$Biorep) %>%
  #           map(~ as.matrix(.x)) %>%
  #           do.call(cbind,.))
  # }))
                                       
mix_params <- get_mixing_params(wide_data)
pSups <- calculate_pSup(wide_data,mix_params)

