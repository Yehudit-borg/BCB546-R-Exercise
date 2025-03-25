#loading essential libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

#DATA INSPECTION 
fang <- read.table("fang_et_al_genotypes.txt", header = TRUE, sep = "\t")

# Read the file
transG1 <- read.table("transposed_genotypes.txt", header = TRUE, sep = "\t")
snp <- read.table("snp_position.txt", header = TRUE, sep = "\t")
dim(snp)
colnames(snp)

str(transG1)
dim(transG1)
colnames(snp)
colnames(transG1)
colnames(fang)
str(fang)
dim(fang)
head(snp)
head(transG1)
head(fang)
tail(snp)
tail(transG1)
tail(fang)

summary(fang)
summary(snp)
summary(transG1)



# Extract the 2nd row as the new header
transNewH <- transG1[2, ]
str(transNewH)
# Remove the original header and the 2nd row
transNewH1 <- transG1[-c(1, 2), ]

# Convert the new header to a character vector
new_header <- as.character(transNewH)

# Assign the new header to the data frame
colnames(transNewH1) <- new_header

# View the updated data frame

dim(transNewH1)


colnames(transNewH1)
# Rename the first column
colnames(transNewH1)[1] <- "SNP_ID"
colnames(transNewH1)
dim(transNewH1)

# Make column names unique
colnames(transNewH1) <- make.unique(colnames(transNewH1))

# Verify the updated column names
colnames(transNewH1)

#Merging the snp and transGSNP datasets
merged_data <- full_join(snp, transNewH1, by = "SNP_ID")
dim(merged_data)
colnames(merged_data)

#Extracting the column  names to identify the regions where each group appears
colname <- colnames(merged_data)
#Saving the group names (column names) in a csv file
write.csv(colname, "output.csv", row.names=FALSE, quote=FALSE) 

#Extracting required columns (Subsetting)
merged_extract <- merged_data %>%
  select(1, 3, 4, 2508:2797, 1225:2480, 2481:2507)

Maize_increase_numeric_positions <- merged_extract[order(as.numeric(merged_extract$Position)), ]

unique(merged_extract$Position)

# Extract and sort for Chromosomes 
extract_and_sort <- function(chromosome_number, data, order_type = "asc") {
  sorted_data <- data %>%
    filter(Chromosome == chromosome_number) %>%
    mutate(Position = as.numeric(Position)) %>%
    arrange(if (order_type == "asc") Position else desc(Position))
  
  return(sorted_data)
}

# Define a vector of chromosome numbers
chromosomes <- 1:10  

# Arranging in ascending order
# Use lapply to extract and sort all chromosomes at once in ascending order
maize_asc_list <- lapply(chromosomes, extract_and_sort, data = merged_extract, order_type = "asc")

# Assign names to the list for easy reference
names(maize_asc_list) <- paste0("maize", chromosomes, "_asc")


# View the extracted and sorted chromosome data
maize_asc_list

# Arranging in ascending order
# Use lapply to extract and sort all chromosomes at once in descendinfg order
maize_desc_list <- lapply(chromosomes, extract_and_sort, data = merged_extract, order_type = "desc")

# Assign names to the list for easy reference
names(maize_desc_list) <- paste0("maize", chromosomes, "_desc")

# View the extracted and sorted chromosome data
maize_desc_list


#Saving the files
# Create directories to store the files
dir.create("maize_asc", showWarnings = FALSE)
dir.create("maize_desc", showWarnings = FALSE)

# Function to save files
save_csv <- function(data, name, folder) {
  file_name <- paste0(folder, "/", name, ".csv")
  write.csv(data, file = file_name, row.names = FALSE)
}

# Save ascending files using lapply
lapply(names(maize_asc_list), function(chr) {
  save_csv(maize_asc_list[[chr]], paste0(chr, "_sortasc"), "maize_asc")
})

# Save descending files using lapply
lapply(names(maize_desc_list), function(chr) {
  save_csv(maize_desc_list[[chr]], paste0(chr, "_sortdesc"), "maize_desc")
})





#Extracting required Teosinte columns (Subsetting)
Teosinte_extract <- merged_data %>%
  select(1, 3, 4, 89:988, 1178:1218, 989:1022)

Teosinte_numeric_positions <- Teosinte_extract[order(as.double(Teosinte_extract$Position)), ]

unique(Teosinte_extract$Position)

# Extract and sort for Chromosomes
Teoextract_and_sort <- function(chromosome_number, data, order_type = "asc") {
  sorted_data <- data %>%
    filter(Chromosome == chromosome_number) %>%
    mutate(Position = as.numeric(Position)) %>%
    arrange(if (order_type == "asc") Position else desc(Position))
  
  return(sorted_data)
}

# Use lapply to extract and sort all chromosomes at once in ascending order
Teo_asc_list <- lapply(chromosomes, Teoextract_and_sort, data = Teosinte_extract, order_type = "asc")

# Assign names to the list for easy reference
names(Teo_asc_list) <- paste0("teosinte", chromosomes, "_asc")


# View the extracted and sorted chromosome data
Teo_asc_list

# Use lapply to extract and sort all chromosomes at once in descending order
Teo_desc_list <- lapply(chromosomes, Teoextract_and_sort, data = Teosinte_extract, order_type = "desc")

# Assign names to the list for easy reference
names(Teo_desc_list) <- paste0("teosinte", chromosomes, "_desc")

# View the extracted and sorted chromosome data
Teo_desc_list

# Create directories to store the files
dir.create("teosinte_asc", showWarnings = FALSE)
dir.create("teosinte_desc", showWarnings = FALSE)

# Function to save files
save_csv <- function(data, name, folder) {
  file_name <- paste0(folder, "/", name, ".csv")
  write.csv(data, file = file_name, row.names = FALSE)
}

# Save ascending files using lapply
lapply(names(Teo_asc_list), function(chr) {
  save_csv(Teo_asc_list[[chr]], paste0(chr, "_sortasc"), "teosinte_asc")
})

# Save descending files using lapply
lapply(names(Teo_desc_list), function(chr) {
  save_csv(Teo_desc_list[[chr]], paste0(chr, "_sortdesc"), "teosinte_desc")
})



# Clean and prepare Maize data
merged_extract_clean <- merged_extract %>%
  drop_na(Chromosome, Position) %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  filter(Chromosome != "" & !is.na(as.numeric(Chromosome))) %>%
  mutate(Chromosome = factor(Chromosome, levels = sort(as.numeric(unique(Chromosome))))) %>%
  pivot_longer(cols = starts_with("SNP_"), names_to = "SNP_ID", values_to = "Genotype") %>%
  drop_na(Genotype) %>%
  group_by(Chromosome) %>%
  summarise(SNP_Count = n(), .groups = "drop") %>%
  mutate(Species = "Maize")  # Label dataset

# Clean and prepare Teosinte data
Teosinte_extract_clean <- Teosinte_extract %>%
  drop_na(Chromosome, Position) %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  filter(Chromosome != "" & !is.na(as.numeric(Chromosome))) %>%
  mutate(Chromosome = factor(Chromosome, levels = sort(as.numeric(unique(Chromosome))))) %>%
  pivot_longer(cols = starts_with("SNP_"), names_to = "SNP_ID", values_to = "Genotype") %>%
  drop_na(Genotype) %>%
  group_by(Chromosome) %>%
  summarise(SNP_Count = n(), .groups = "drop") %>%
  mutate(Species = "Teosinte")  # Label dataset

# Combine both datasets
combined_snp_counts <- bind_rows(merged_extract_clean, Teosinte_extract_clean)

# Composite bar plot comparing SNP counts for Maize and Teosinte
composite_snp_plot <- ggplot(combined_snp_counts, aes(x = Chromosome, y = SNP_Count, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +  # "Dodge" for side-by-side bars
  theme_minimal() +
  labs(title = "Comparison of SNP Counts per Chromosome: Maize vs Teosinte",
       x = "Chromosome",
       y = "Number of SNPs",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Display the plot
print(composite_snp_plot)

# Save the plot as a PDF
ggsave("Maize_vs_Teosinte_SNP_count_per_chromosome.pdf", plot = composite_snp_plot, width = 10, height = 6)





################ DISTRIBUTION OF SNPS ACROSS CHROMOSOME

# Scatter plot of SNPs distribution across chromosomes for Maize
merged_extract_clean1 <- merged_extract %>%
  mutate(Position = as.numeric(Position),  # Ensure Position is numeric
         Chromosome = as.numeric(as.character(Chromosome))) %>%  # Convert Chromosome to numeric
  drop_na(Chromosome, Position) %>%  # Remove missing values
  mutate(Chromosome = factor(Chromosome, levels = sort(unique(Chromosome))))  # Reorder Chromosome levels

dist_plot <- ggplot(merged_extract_clean1, aes(x = Position, y = Chromosome, color = Chromosome)) +
  geom_point(alpha = 0.5, size = 1) +  # Plot points with transparency (alpha = 0.5)
  theme_minimal() +  # Apply minimal theme
  labs(title = "Distribution of SNPs Across Chromosomes (Maize)",  # Add plot title
       x = "Genomic Position",  # Label for the x-axis
       y = "Chromosome") +  # Label for the y-axis
  theme(legend.position = "none")  # Remove legend for cleaner plot
  
(dist_plot)

# Save the SNP distribution plot as a PDF
ggsave("Maize_SNP_distribution_across_chromosomes.pdf", plot = dist_plot, width = 8, height = 6)


# Scatter plot of SNPs distribution across chromosomes for Teosinte
Teosinte_extract_clean1 <- Teosinte_extract %>%
  mutate(Position = as.numeric(Position),  # Ensure Position is numeric
         Chromosome = as.numeric(as.character(Chromosome))) %>%  # Convert Chromosome to numeric
  drop_na(Chromosome, Position) %>%  # Remove missing values
  mutate(Chromosome = factor(Chromosome, levels = sort(unique(Chromosome))))  # Reorder Chromosome levels

Teo_dist_plot <- ggplot(Teosinte_extract_clean1, aes(x = Position, y = Chromosome, color = Chromosome)) +
  geom_point(alpha = 0.5, size = 1) +  # Plot points with transparency (alpha = 0.5)
  theme_minimal() +  # Apply minimal theme
  labs(title = "Distribution of SNPs Across Chromosomes (Teosinte)",  # Add plot title
       x = "Genomic Position",  # Label for the x-axis
       y = "Chromosome") +  # Label for the y-axis
  theme(legend.position = "none")  # Remove legend for cleaner plot
(Teo_dist_plot)

# Save the SNP distribution plot as a PDF
ggsave("Teosinte_SNP_distribution_across_chromosomes.pdf", plot = Teo_dist_plot, width = 8, height = 6)


###### SAMPLES: HOMOGENEOUS-HETEROGENOUS-MISSING DATA PROPORTION

# Reshape the data to long format
fang_long <- fang %>%
  pivot_longer(cols = -c(Sample_ID, JG_OTU, Group), names_to = "Site", values_to = "Genotype")

# Function to classify genotypes
classify_genotype <- function(geno) {
  if (geno == "?/?") {
    return("Missing")
  } else if (substr(geno, 1, 1) == substr(geno, 3, 3)) {
    return("Homozygous")
  } else {
    return("Heterozygous")
  }
}

# Apply classification
fang_long <- fang_long %>%
  mutate(Genotype_Type = sapply(Genotype, classify_genotype))

# Summarize proportions per Sample_ID
summary_table <- fang_long %>%
  group_by(Sample_ID, Genotype_Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Proportion = Count / sum(Count))  # Calculate proportion directly

# Stacked bar plot (normalized height using "position = fill")
sample_plot <- ggplot(summary_table, aes(x = Sample_ID, y = Proportion, fill = Genotype_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Normalize heights automatically
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  theme_minimal() +
  labs(title = "Proportion of Homozygous, Heterozygous, and Missing Data Per Sample",
       x = "Sample ID",
       y = "Proportion",
       fill = "Genotype Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

# Display the plot
(sample_plot)

# Save the plot
ggsave("Sample_HHM_Proportion.pdf", plot = sample_plot, width = 10, height = 6)



###### GROUPS: HOMOGENEOUS-HETEROGENOUS-MISSING DATA PROPORTION

library(dplyr)
library(tidyr)
library(readr)  # For writing CSV files

# Reshape the data to long format
fang_long <- fang %>%
  pivot_longer(cols = -c(Sample_ID, JG_OTU, Group), names_to = "Site", values_to = "Genotype")

# Function to classify genotypes
classify_genotype <- function(geno) {
  if (geno == "?/?") {
    return("Missing")
  } else if (substr(geno, 1, 1) == substr(geno, 3, 3)) {
    return("Homozygous")
  } else {
    return("Heterozygous")
  }
}

# Apply classification
fang_long <- fang_long %>%
  mutate(Genotype_Type = sapply(Genotype, classify_genotype))

# Summarize counts and proportions per Group
group_summary <- fang_long %>%
  group_by(Group, Genotype_Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Proportion = (Count / sum(Count)) * 100)  # Calculate proportion directly

# View the output table
print(group_summary)

# Convert data to long format for ggplot (directly use the summary without pivot_wider)
group_long <- group_summary %>%
  select(Group, Genotype_Type, Proportion)  # No need for pivot_wider

# Stacked bar plot for groups
group_plot <- ggplot(group_long, aes(x = Group, y = Proportion, fill = Genotype_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Normalize bar height
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  theme_minimal() +
  labs(title = "Proportion of Homozygous, Heterozygous, and Missing Data Per Group",
       x = "Group",
       y = "Proportion",
       fill = "Genotype Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# Display the plot
print(group_plot)

# Save the plot
ggsave("Group_HHM_Proportion.pdf", plot = group_plot, width = 10, height = 6)





# Count the number of Sample_IDs per Group: TABLE
sample_count_per_group <- fang_long %>%
  group_by(Group) %>%
  summarise(Sample_Count = n_distinct(Sample_ID), .groups = "drop")

# Print the result
print(sample_count_per_group)


# Count the number of Sample_IDs per Group: VISUALIZATION
sample_count_per_group <- fang_long %>%
  group_by(Group) %>%
  summarise(Sample_Count = n_distinct(Sample_ID), .groups = "drop")

# Plotting the bar chart
library(ggplot2)
plot <- ggplot(sample_count_per_group, aes(x = Group, y = Sample_Count)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Sample Count per Group", x = "Group", y = "Sample Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
(plot)

# Save the plot as a PDF
ggsave("sample_count_per_group.pdf", plot = plot, width = 8, height = 6)
