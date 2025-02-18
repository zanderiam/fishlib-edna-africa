#!/usr/bin/env Rscript

# R script to make reference databases for UK fishes for multiple markers
source(here::here("scripts/load-libs.R"))
source(here::here("scripts/ncbi-utils.R"))
source(here::here("scripts/data-processing.R"))
source(here::here("scripts/taxonomy-utils.R"))

# Parse command line arguments
option_list <- list(
  make_option(c("-t","--threads"), type="numeric"),
  make_option(c("-m","--metabarcode"), type="character")
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
cores <- min(opt$threads, parallel::detectCores() - 1)

# Load and process initial data
data <- load_initial_data()
prefixes <- process_metabarcode_prefixes(opt$metabarcode)

# Extract fragments using HMMER
writeLines("\nExtracting metabarcode fragments with HMMER (may take several minutes) ...")
fragments <- extract_fragments(prefixes, cores)
writeLines("\nDone")

# Process NCBI data
writeLines("\nRetrieving metadata from NCBI ...")
start_time <- Sys.time()
ncbi_results <- process_ncbi_data(fragments, cores, data$bold_data)
end_time <- Sys.time()
cat("\nTime elapsed:", difftime(end_time, start_time, units = "secs"), "\n")

# Process taxonomy and create final dataset
final_data <- process_taxonomy(ncbi_results, data$species_table, data$bold_data, data$stats)

# Write output
writeLines("\nWriting out reference library to 'assets/reference-library-master.csv.gz' ...")
write_csv(final_data, file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
writeLines("\nAll operations completed!")