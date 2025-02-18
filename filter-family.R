#!/usr/bin/env Rscript

# Load necessary library
library(dplyr)

# Read the CSV file
species_data <- read.csv("assets/species-table.csv")

# Keep one row for each unique family name
unique_families <- species_data %>%
  group_by(family) %>%
  slice(1) %>%
  ungroup()

# Write the result to a new CSV file
write.csv(unique_families, "unique-families.csv", row.names = FALSE)

# Print the result to the console
print(unique_families)
