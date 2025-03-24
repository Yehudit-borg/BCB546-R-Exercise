#loading essential libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

#DATA INSPECTION 
fang <- read.table("fang_et_al_genotypes.txt", header = TRUE, sep = "\t")
View(fang)
# Read the file
transG1 <- read.table("transposed_genotypes.txt", header = TRUE, sep = "\t")
snp <- read.table("snp_position.txt", header = TRUE, sep = "\t")
dim(snp)
colnames(snp)
View(snp)
View(transG1)
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
View(transNewH1)
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
View(merged_data)

#Extracting the column  names to identify the regions where each group appears
colname <- colnames(merged_data)
#Saving the group names (column names) in a csv file
write.csv(colname, "output.csv", row.names=FALSE, quote=FALSE) 

#Extracting required columns (Subsetting)
merged_extract <- merged_data[, c(1, 3, 4, 2508:2797, 1225:2480, 2481:2507)]
View(merged_extract)

#Subsetting by sorting the Maize file in increasing and decreasing order
Maize_descrease <- arrange(merged_extract, desc(Position))
View(Maize_descrease)

Maize_increase <- arrange(merged_extract, Position)
View(Maize_increase)

Maize_increase_numeric_positions <- merged_extract[order(as.double(merged_extract$Position)), ]
View(Maize_increase_numeric_positions)
unique(merged_extract$Position)


# Function to extract chromosome and sort by Position
extract_and_sort <- function(chromosome_number, data) {
  data %>%
    filter(.[[2]] == chromosome_number) %>%  # Extract specific chromosome
    mutate(Position_numeric = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric
    arrange(Position_numeric) %>%  # Sort in increasing order
    select(-Position_numeric)  # Remove the helper column
}

# Extract and sort for Chromosome 1
chr1_asc <- extract_and_sort(1, merged_extract)

chr1_positions <- chr1_asc$Position  # Extract Position column
head(chr1_positions)  # Check first few values
tail(chr1_positions)  # Check last few values

print(chr1_positions)

is_sorted <- all(diff(as.numeric(chr1_positions)) >= 0)
print(is_sorted) 

# Extract and sort for Chromosome 2
chr2_asc <- extract_and_sort(2, merged_extract)

# Extract and sort for Chromosome 3
chr3_asc <- extract_and_sort(3, merged_extract)

# Extract and sort for Chromosome 4
chr4_asc <- extract_and_sort(4, merged_extract)

# Extract and sort for Chromosome 5
chr5_asc <- extract_and_sort(5, merged_extract)

# Extract and sort for Chromosome 6
chr6_asc <- extract_and_sort(6, merged_extract)

# Extract and sort for Chromosome 7
chr7_asc <- extract_and_sort(7, merged_extract)

# Extract and sort for Chromosome 8
chr8_asc <- extract_and_sort(8, merged_extract)

# Extract and sort for Chromosome 9
chr9_asc <- extract_and_sort(9, merged_extract)

# Extract and sort for Chromosome 10
chr10_asc <- extract_and_sort(10, merged_extract)

extract_and_sort_desc <- function(chromosome_number, data) {
  data %>%
    filter(.[[2]] == chromosome_number) %>%  # Extract specific chromosome
    mutate(Position_numeric = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric
    arrange(desc(Position_numeric)) %>%  # Sort in descending order
    mutate(across(where(is.character), ~ gsub("\\?", "-", .))) %>%  # Replace "?" with "-" in character columns
    select(-Position_numeric)  # Drop the helper column
}

chr1_desc <- extract_and_sort_desc(1, merged_extract)
chr2_desc <- extract_and_sort_desc(2, merged_extract)
chr3_desc <- extract_and_sort_desc(3, merged_extract)
chr4_desc <- extract_and_sort_desc(4, merged_extract)
chr5_desc <- extract_and_sort_desc(5, merged_extract)
chr6_desc <- extract_and_sort_desc(6, merged_extract)
chr7_desc <- extract_and_sort_desc(7, merged_extract)
chr8_desc <- extract_and_sort_desc(8, merged_extract)
chr9_desc <- extract_and_sort_desc(9, merged_extract)
chr10_desc <- extract_and_sort_desc(10, merged_extract)


# Create a directory to store the files in ascending order
dir.create("chromosome_asc", showWarnings = FALSE)

# List of chromosome datasets
chromosome_asclist <- list(
  chr1 = chr1_asc, chr2 = chr2_asc, chr3 = chr3_asc,
  chr4 = chr4_asc, chr5 = chr5_asc, chr6 = chr6_asc,
  chr7 = chr7_asc, chr8 = chr8_asc, chr9 = chr9_asc,
  chr10 = chr10_asc
)

# Loop through and save each chromosome dataset
for (chr in names(chromosome_asclist)) {
  file_name <- paste0("chromosome_asc/", chr, "_sortasc.csv")
  write.csv(chromosome_asclist[[chr]], file = file_name, row.names = FALSE)
}



# Create a directory to store the files in descending order
dir.create("chromosome_desc", showWarnings = FALSE)

# List of chromosome datasets
chromosome_desclist <- list(
  chr1 = chr1_desc, chr2 = chr2_desc, chr3 = chr3_desc,
  chr4 = chr4_desc, chr5 = chr5_desc, chr6 = chr6_desc,
  chr7 = chr7_desc, chr8 = chr8_desc, chr9 = chr9_desc,
  chr10 = chr10_desc
)

# Loop through and save each chromosome dataset
for (chr in names(chromosome_desclist)) {
  file_name <- paste0("chromosome_desc/", chr, "_sortdesc.csv")
  write.csv(chromosome_desclist[[chr]], file = file_name, row.names = FALSE)
}





#Extracting the column  names to identify the regions where each group appears
colname <- colnames(merged_data)
#Saving the group names (column names) in a csv file
write.csv(colname, "output.csv", row.names=FALSE, quote=FALSE) 

#Extracting required Teosinte columns (Subsetting)
Teosinte_extract <- merged_data[, c(1, 3, 4, 89:988, 1178:1218, 989:1022)]
View(Teosinte_extract)


Teosinte_numeric_positions <- Teosinte_extract[order(as.double(Teosinte_extract$Position)), ]
View(Teosinte_numeric_positions)
unique(Teosinte_extract$Position)


# Function to extract chromosome and sort by Position
Teoextract_and_sort <- function(chromosome_number, data) {
  data %>%
    filter(.[[2]] == chromosome_number) %>%  # Extract specific chromosome
    mutate(Position_numeric = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric
    arrange(Position_numeric) %>%  # Sort in increasing order
    select(-Position_numeric)  # Remove the helper column
}

# Extract and sort for Chromosome 1
Tchr1_asc <- Teoextract_and_sort(1, Teosinte_extract)

Tchr1_positions <- chr1_asc$Position  # Extract Position column
head(Tchr1_positions)  # Check first few values
tail(Tchr1_positions)  # Check last few values

print(chr1_positions)

Tis_sorted <- all(diff(as.numeric(Tchr1_positions)) >= 0)
print(Tis_sorted) 

# Extract and sort for Chromosome 2
Tchr2_asc <- Teoextract_and_sort(2, Teosinte_extract)

# Extract and sort for Chromosome 3
Tchr3_asc <- Teoextract_and_sort(3, Teosinte_extract)

# Extract and sort for Chromosome 4
Tchr4_asc <- Teoextract_and_sort(4, Teosinte_extract)

# Extract and sort for Chromosome 5
Tchr5_asc <- Teoextract_and_sort(5, Teosinte_extract)

# Extract and sort for Chromosome 6
Tchr6_asc <- Teoextract_and_sort(6, Teosinte_extract)

# Extract and sort for Chromosome 7
Tchr7_asc <- Teoextract_and_sort(7, Teosinte_extract)

# Extract and sort for Chromosome 8
Tchr8_asc <- Teoextract_and_sort(8, Teosinte_extract)

# Extract and sort for Chromosome 9
Tchr9_asc <- Teoextract_and_sort(9, Teosinte_extract)

# Extract and sort for Chromosome 10
Tchr10_asc <- Teoextract_and_sort(10, Teosinte_extract)

Teoextract_and_sort_desc <- function(chromosome_number, data) {
  data %>%
    filter(.[[2]] == chromosome_number) %>%  # Extract specific chromosome
    mutate(Position_numeric = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric
    arrange(desc(Position_numeric)) %>%  # Sort in descending order
    mutate(across(where(is.character), ~ gsub("\\?", "-", .))) %>%  # Replace "?" with "-" in character columns
    select(-Position_numeric)  # Drop the helper column
}

Tchr1_desc <- Teoextract_and_sort_desc(1, merged_extract)
Tchr2_desc <- Teoextract_and_sort_desc(2, merged_extract)
Tchr3_desc <- Teoextract_and_sort_desc(3, merged_extract)
Tchr4_desc <- Teoextract_and_sort_desc(4, merged_extract)
Tchr5_desc <- Teoextract_and_sort_desc(5, merged_extract)
Tchr6_desc <- Teoextract_and_sort_desc(6, merged_extract)
Tchr7_desc <- Teoextract_and_sort_desc(7, merged_extract)
Tchr8_desc <- Teoextract_and_sort_desc(8, merged_extract)
Tchr9_desc <- Teoextract_and_sort_desc(9, merged_extract)
Tchr10_desc <- Teoextract_and_sort_desc(10, merged_extract)


# Create a directory to store the files in ascending order
dir.create("Teochromosome_asc", showWarnings = FALSE)

# List of chromosome datasets
Tchromosome_asclist <- list(
  chr1 = Tchr1_asc, chr2 = Tchr2_asc, chr3 = Tchr3_asc,
  chr4 = Tchr4_asc, chr5 = Tchr5_asc, chr6 = Tchr6_asc,
  chr7 = Tchr7_asc, chr8 = Tchr8_asc, chr9 = Tchr9_asc,
  chr10 = Tchr10_asc
)

# Loop through and save each chromosome dataset
for (chr in names(Tchromosome_asclist)) {
  file_name <- paste0("Teochromosome_asc/", chr, "_teoasc.csv")
  write.csv(Tchromosome_asclist[[chr]], file = file_name, row.names = FALSE)
}


# Create a directory to store the files in descending order
dir.create("Teochromosome_desc", showWarnings = FALSE)

# List of chromosome datasets
Tchromosome_desclist <- list(
  chr1 = Tchr1_desc, chr2 = Tchr2_desc, chr3 = Tchr3_desc,
  chr4 = Tchr4_desc, chr5 = Tchr5_desc, chr6 = Tchr6_desc,
  chr7 = Tchr7_desc, chr8 = Tchr8_desc, chr9 = Tchr9_desc,
  chr10 = Tchr10_desc
)

# Loop through and save each chromosome dataset
for (chr in names(Tchromosome_desclist)) {
  file_name <- paste0("Teochromosome_desc/", chr, "_teodesc.csv")
  write.csv(Tchromosome_desclist[[chr]], file = file_name, row.names = FALSE)
}



# Count SNPs per chromosome for Maize
# Remove rows with missing values in Chromosome or Position
merged_extract_clean <- merged_extract %>%
  drop_na(Chromosome, Position)

# Check for NA values in the Chromosome column
sum(is.na(merged_extract_clean$Chromosome))


# Check unique values in the Chromosome column
unique(merged_extract_clean$Chromosome)

# Remove rows with NA, blank, or non-numeric Chromosome values
merged_extract_clean <- merged_extract %>%
  filter(!is.na(Chromosome) & Chromosome != "" & !is.na(as.numeric(as.character(Chromosome))))

# Check for NA values in the Chromosome column
sum(is.na(merged_extract_clean$Chromosome))


# Ensure Chromosome is treated as a numeric factor in sequential order
merged_extract_clean$Chromosome <- factor(
  merged_extract_clean$Chromosome,
  levels = sort(as.numeric(unique(merged_extract_clean$Chromosome)))
)

# Count SNPs per chromosome
snp_counts <- merged_extract_clean %>%
  group_by(Chromosome) %>%  # Group data by chromosome
  summarise(SNP_Count = n())  # Count the number of SNPs for each chromosome

# Bar plot of SNP counts per chromosome
snp_plot <- ggplot(snp_counts, aes(x = as.factor(Chromosome), y = SNP_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Create bar plot with steelblue color
  theme_minimal() +  # Apply minimal theme for cleaner visuals
  labs(title = "Number of SNPs per Chromosome (Maize)",  # Add plot title
       x = "Chromosome",  # Label for the x-axis
       y = "Number of SNPs") +  # Label for the y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

(snp_plot)

# Save the SNP count plot as a PDF
ggsave("Maize_SNP_count_per_chromosome.pdf", plot = snp_plot, width = 8, height = 6)



# Count SNPs per chromosome for Teosinte

# Remove rows with missing values in Chromosome or Position
Teosinte_extract_clean <- Teosinte_extract %>%
  drop_na(Chromosome, Position)

# Check for NA values in the Chromosome column
sum(is.na(Teosinte_extract_clean$Chromosome))


# Check unique values in the Chromosome column
unique(Teosinte_extract_clean$Chromosome)

# Remove rows with NA, blank, or non-numeric Chromosome values
Teosinte_extract_clean <- Teosinte_extract %>%
  filter(!is.na(Chromosome) & Chromosome != "" & !is.na(as.numeric(as.character(Chromosome))))

# Check for NA values in the Chromosome column
sum(is.na(Teosinte_extract_clean$Chromosome))


# Ensure Chromosome is treated as a numeric factor in sequential order
Teosinte_extract_clean$Chromosome <- factor(
  Teosinte_extract_clean$Chromosome,
  levels = sort(as.numeric(unique(Teosinte_extract_clean$Chromosome)))
)

# Count SNPs per chromosome
Teo_snp_counts <- Teosinte_extract_clean %>%
  group_by(Chromosome) %>%  # Group data by chromosome
  summarise(SNP_Count = n())  # Count the number of SNPs for each chromosome

# Bar plot of SNP counts per chromosome
Teo_snp_plot <- ggplot(Teo_snp_counts, aes(x = Chromosome, y = SNP_Count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +  # Create bar plot with steelblue color
  theme_minimal() +  # Apply minimal theme for cleaner visuals
  labs(title = "Number of SNPs per Chromosome (Teosinte)",  # Add plot title
       x = "Chromosome",  # Label for the x-axis
       y = "Number of SNPs") +  # Label for the y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Display the plot
print(Teo_snp_plot)


# Save the SNP count plot as a PDF
ggsave("Teosinte_SNP_count_per_chromosome.pdf", plot = Teo_snp_plot, width = 8, height = 6)




# # Add species column to each dataset
 maize_snp_counts <- merged_extract_clean %>%
   group_by(Chromosome) %>%
   summarise(SNP_Count = n()) %>%
   mutate(Species = "Maize")  # Label as Maize
 
 teosinte_snp_counts <- Teosinte_extract_clean %>%
   group_by(Chromosome) %>%
   summarise(SNP_Count = n()) %>%
   mutate(Species = "Teosinte")  # Label as Teosinte

# Combine both datasets into one
combined_snp_counts <- bind_rows(maize_snp_counts, teosinte_snp_counts)

# Plot SNP counts per chromosome for Maize and Teosinte
comparison_plot <- ggplot(combined_snp_counts, aes(x = as.factor(Chromosome), y = SNP_Count, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge for side-by-side bars
  theme_minimal() +  # Clean theme
  scale_fill_manual(values = c("Maize" = "steelblue", "Teosinte" = "darkgreen")) +  # Custom colors
  labs(title = "Comparison of SNPs per Chromosome (Maize vs. Teosinte)",
       x = "Chromosome",
       y = "Number of SNPs",
       fill = "Species") +  # Legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# Display the plot
(comparison_plot)

# Save the plot
ggsave("Maize_vs_Teosinte_SNP_count_per_chromosome.pdf", plot = comparison_plot, width = 8, height = 6)




################ DISTRIBUTION OF SNPS ACROSS CHROMOSOME

# Scatter plot of SNPs distribution across chromosomes for Maize
dist_plot <- ggplot(merged_extract_clean, aes(x = Position, y = Chromosome, color = Chromosome)) +
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
Teo_dist_plot <- ggplot(Teosinte_extract_clean, aes(x = Position, y = Chromosome, color = Chromosome)) +
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
  group_by(Sample_ID) %>%
  mutate(Proportion = (Count / sum(Count)) * 100) %>%  # Convert to percentage
  select(Sample_ID, Genotype_Type, Proportion) %>%
  pivot_wider(names_from = Genotype_Type, values_from = Proportion, values_fill = 0) %>%  # Fill missing with 0
  mutate(across(where(is.numeric), round, 2))  # Round to two decimal places


# View the output table
print(summary_table)


# Convert data to long format for ggplot
sample_long <- summary_table %>%
  pivot_longer(cols = -Sample_ID, names_to = "Genotype_Type", values_to = "Proportion")

# Stacked bar plot
sample_plot <- ggplot(sample_long, aes(x = Sample_ID, y = Proportion, fill = Genotype_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot with normalized height
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

# Summarize proportions per Group
group_summary <- fang_long %>%
  group_by(Group, Genotype_Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Proportion = (Count / sum(Count)) * 100) %>%  # Convert to percentage
  select(Group, Genotype_Type, Proportion) %>%
  pivot_wider(names_from = Genotype_Type, values_from = Proportion, values_fill = 0) %>%  # Fill missing with 0
  mutate(across(where(is.numeric), round, 2))  # Round to two decimal places


# View the output table
print(group_summary)


# Convert data to long format for ggplot
group_long <- group_summary %>%
  pivot_longer(cols = -Group, names_to = "Genotype_Type", values_to = "Proportion")

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
(group_plot)

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
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Sample Count per Group", x = "Group", y = "Sample Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Save the plot as a PDF
ggsave("sample_count_per_group.pdf", plot = plot, width = 8, height = 6)
