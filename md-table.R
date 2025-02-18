#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read and process the markdown table
data <- read.table("reports/reports-table-only.md", 
                   header = FALSE, 
                   sep = "|", 
                   skip = 0,
                   strip.white = TRUE,
                   check.names = FALSE) %>%
  # Remove first and last empty columns from pipe table
  select(-1, -ncol(.))

# Set column names from the first row and remove it
colnames(data) <- as.character(unlist(data[1,]))
data <- data[-1,] # Remove the header row
data <- data[-1,] # Remove the separator row

# Clean column names
colnames(data) <- c("Locus", "Fragment_Name", "Fragment_Size", "Total", "Cov_all")

# Convert numeric columns, keeping Fragment_Name as character
data$Fragment_Size <- as.numeric(gsub("[^0-9.]", "", data$Fragment_Size))
data$Total <- as.numeric(gsub("[^0-9.]", "", data$Total))
data$Cov_all <- as.numeric(gsub("[^0-9.]", "", data$Cov_all))

# Create a directory for plots if it doesn't exist
dir.create("plots", showWarnings = FALSE)

# Create separate plots for each locus
# unique_loci <- unique(data$Locus)
# 
# for (locus in unique_loci) {
#   locus_data <- data %>% filter(Locus == locus)
  
  plot <- ggplot(data, aes(x = Fragment_Name, y = Total)) +
    geom_bar(stat = "identity", fill = "#2E86C1") +
    labs(
      title = paste("Coverage Analysis -", Fragment_Name),
      subtitle = "Total Count per Fragment",
      x = "Fragment",
      y = "Coverage"
    ) +
    theme_economist_white() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "#E0E0E0"),
      panel.grid.minor = element_line(color = "#F5F5F5")
    )
  
  # Save the plot with better resolution
  ggsave(
    filename = file.path("plots", paste0(Fragment_Name, "_coverage_analysis.png")),
    plot = plot,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  cat(sprintf("Plot saved: %s_coverage_analysis.png\n",Fragment_Name))
#}

# Print summary statistics
cat("\nSummary Statistics:\n")
summary_stats <- data %>%
  group_by(Locus) %>%
  summarise(
    Mean_Coverage = mean(Total),
    Min_Coverage = min(Total),
    Max_Coverage = max(Total),
    n_Fragments = n()
  )
print(summary_stats)
