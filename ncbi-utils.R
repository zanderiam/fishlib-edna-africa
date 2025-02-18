#!/usr/bin/env Rscript

# Configuration
MAX_RETRIES <- 5          # Increased retries
RETRY_DELAY <- 10         # Increased delay between retries
QUERY_DELAY <- 1          # Increased delay between queries
CHUNK_SIZE <- 100         # Reduced chunk size for better reliability

#' Safe NCBI query with retries and detailed error reporting
#' @param chunk Vector of accession numbers
#' @return Data frame of results or NULL on failure
ncbi_query_with_retry <- function(chunk) {
  for (attempt in 1:MAX_RETRIES) {
    tryCatch({
      # Add delay between retries
      if (attempt > 1) {
        cat(sprintf("\nRetry attempt %d/%d after %d seconds...", attempt, MAX_RETRIES, RETRY_DELAY))
        Sys.sleep(RETRY_DELAY)
      }
      
      result <- rentrez::entrez_fetch(
        db = "nucleotide",
        id = chunk,
        rettype = "gb",
        retmode = "xml"
      )
      
      if (is.null(result) || length(result) == 0) {
        warning("Empty result received from NCBI")
        next
      }
      
      # Parse and validate result
      parsed <- parse_ncbi_result(result)
      if (is.null(parsed) || nrow(parsed) == 0) {
        warning("Failed to parse NCBI result")
        next
      }
      
      return(parsed)
      
    }, error = function(e) {
      warning(sprintf("NCBI query failed: %s", e$message))
      return(NULL)
    })
  }
  return(NULL)
}

#' Parse NCBI result with error handling
#' @param result Raw XML result from NCBI
#' @return Parsed data frame or NULL on failure
parse_ncbi_result <- function(result) {
  tryCatch({
    # Parse XML and convert to data frame
    # Add your existing parsing logic here
    parsed <- xml2::read_xml(result)
    # ... rest of parsing logic
    return(data.frame(...))
  }, error = function(e) {
    warning(sprintf("Failed to parse NCBI result: %s", e$message))
    return(NULL)
  })
}

#' Process NCBI queries in batches with progress tracking and error recovery
#' @param accessions Vector of accession numbers
#' @param chunk_size Size of each chunk
#' @param cores Number of cores to use
#' @return List of data frames with results
process_ncbi_queries <- function(accessions, chunk_size = CHUNK_SIZE, cores = 1) {
  # Split into smaller chunks
  chunks <- split(accessions, ceiling(seq_along(accessions)/chunk_size))
  total_chunks <- length(chunks)
  successful_chunks <- 0
  failed_chunks <- 0
  
  # Process chunks with progress reporting
  results <- mcmapply(
    FUN = function(chunk) {
      chunk_num <- match(chunk[1], unlist(chunks))
      cat(sprintf("\rProcessing chunk %d/%d (Success: %d, Failed: %d)...", 
                  chunk_num, total_chunks, successful_chunks, failed_chunks))
      
      Sys.sleep(QUERY_DELAY)
      
      result <- ncbi_query_with_retry(chunk)
      
      if (is.null(result)) {
        failed_chunks <<- failed_chunks + 1
        warning(sprintf("\nChunk %d failed after all retries", chunk_num))
      } else {
        successful_chunks <<- successful_chunks + 1
      }
      
      return(result)
    },
    chunks,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE,
    mc.cores = cores
  )
  
  # Final status report
  cat(sprintf("\nCompleted: %d successful, %d failed chunks\n", 
              successful_chunks, failed_chunks))
  
  # Filter out NULL results
  results <- Filter(Negate(is.null), results)
  
  if (length(results) == 0) {
    stop("All NCBI queries failed. Please check your internet connection and try again.")
  }
  
  if (failed_chunks > 0) {
    warning(sprintf("%d chunks failed. Results may be incomplete.", failed_chunks))
  }
  
  return(results)
}