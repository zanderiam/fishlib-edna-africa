#!/usr/bin/env Rscript
# Load necessary library
library(dplyr)

# Read the dataset
species_data <- read.csv("cleaned-species-table.csv")

# Filter out rows with duplicate speciesName
filtered_data <- species_data %>%
  distinct(fbSpecCode, .keep_all = TRUE)

# Save the filtered dataset to a new CSV file
write.csv(filtered_data, "filtered-species-table.csv", row.names = FALSE)

# Print a message to confirm the operation
cat("Filtered dataset saved as 'filtered-species-table.csv'\n")