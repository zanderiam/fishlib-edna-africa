#!/usr/bin/env Rscript

load_initial_data <- function() {
  tryCatch({
    list(
      species_table = read_csv(file=here("assets/species-table.csv"), show_col_types=FALSE),
      bold_data = read_csv(file=here("temp/bold-dump.csv"), guess_max=100000, show_col_types=FALSE),
      stats = read_csv(file=here("reports/stats.csv"), show_col_types=FALSE)
    )
  }, error = function(e) {
    stop(sprintf("Failed to load initial data: %s", e$message))
  })
}

process_metabarcode_prefixes <- function(metabarcode) {
  prefixes_list <- c("coi.lerayxt","coi.ward","12s.miya","12s.riaz","12s.valentini",
                     "12s.taberlet","16s.berry","cytb.minamoto","16s.kitano","12S.12S-V5", "12s.kelly")
  prefixes_chosen <- unlist(str_split(metabarcode,","))
  
  if(metabarcode == "all") {
    return(paste0(prefixes_list, ".noprimers"))
  } else if (all(prefixes_chosen %in% prefixes_list)) {
    return(paste(prefixes_chosen, "noprimers", sep="."))
  } else {
    stop("'-m' value must be metabarcode(s) listed in Table 1, and separated by a comma")
  }
}

extract_fragments <- function(prefixes, cores) {
  tryCatch({
    fragments <- lapply(prefixes, function(x) {
      run_hmmer3(dir="temp", infile="mtdna-dump.fas", prefix=x, evalue="10", coords="env")
    })
    
    fragments_cat <- do.call(c, fragments)
    fragment_names <- unique(labels(fragments_cat))
    
    if (length(fragment_names) == 0) {
      stop("No fragments extracted")
    }
    
    return(list(
      data = fragments_cat,
      names = fragment_names
    ))
  }, error = function(e) {
    stop(sprintf("Fragment extraction failed: %s", e$message))
  })
}

process_ncbi_data <- function(fragments, cores, bold_data) {
  # Separate GenBank and BOLD data
  in_bold <- fragments$names[fragments$names %in% bold_data$processidUniq]
  in_gb <- fragments$names[!fragments$names %in% bold_data$processidUniq]
  
  if (length(in_gb) == 0 && length(in_bold) == 0) {
    stop("No valid accessions found in either GenBank or BOLD")
  }
  
  # Clean and process GenBank accessions
  in_gb_clean <- in_gb[nchar(in_gb) <= 11]
  if (length(in_gb_clean) == 0) {
    warning("No valid GenBank accessions found")
    return(list(
      ncbi_data = data.frame(),
      bold_ids = in_bold,
      genbank_ids = character(0)
    ))
  }
  
  set.seed(42)
  in_gb_sample <- sample(in_gb_clean)
  
  # Query NCBI with error handling
  tryCatch({
    ncbi_results <- process_ncbi_queries(in_gb_sample, cores=cores)
    if (length(ncbi_results) > 0) {
      frag_df <- as_tibble(bind_rows(ncbi_results))
      writeLines("\nNCBI metadata successfully retrieved")
      return(list(
        ncbi_data = frag_df,
        bold_ids = in_bold,
        genbank_ids = in_gb_clean
      ))
    } else {
      stop("No valid results returned from NCBI queries")
    }
  }, error = function(e) {
    stop(paste("Error processing NCBI queries:", e$message))
  })
}