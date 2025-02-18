#!/usr/bin/env Rscript
# script to quality control the reference libraries and identify erroneous sequences.
# Modified to handle BOLD data separately

# start timer 
start_time <- Sys.time()

# load functions and libs
source(here::here("scripts","load-libs.R"))
source(here::here("scripts","references-load-local.R"))
source(here::here("scripts","references-clean.R"))

# get args
option_list <- list( 
  make_option(c("-t","--threads"), type="numeric"),
  make_option(c("-v","--verbose"), type="character")
)

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

#opt <- NULL
#opt$threads <- 1
#opt$verbose <- "true"

# load stats
stats <- suppressMessages(read_csv(here("reports","stats.csv")))

# message 
writeLines("\nGenerating phylogenetic trees for BOLD data, may take several hours ...")

# set cores
cores <- opt$threads

# make a copy so don't have to reload orig
reflib <- reflib.cleaned

# Separate BOLD and GenBank data
reflib.bold <- reflib %>% filter(source == "BOLD")

# get the prefixes 
prefixes <- reflib.bold %>% 
  select(starts_with("nucleotidesFrag")) %>% 
  names() %>% 
  str_replace_all("nucleotidesFrag\\.","")

# subset each marker for BOLD data
reflibs.sub <- mcmapply(
  function(x) subset_nucs(pref=x, df=reflib.bold), 
  prefixes, 
  SIMPLIFY=FALSE,
  USE.NAMES=TRUE,
  mc.cores=cores
)

# collapse dataframe by haps-per-species, annotate with number haps
reflibs.haps <- mcmapply(
  function(x) haps2fas(df=x), 
  reflibs.sub, 
  SIMPLIFY=FALSE,
  USE.NAMES=TRUE,
  mc.cores=cores
)

# convert to fasta
reflibs.fas <- mcmapply(
  function(x) tab2fas(df=x, seqcol="nucleotidesFrag", namecol="noms"), 
  reflibs.haps, 
  SIMPLIFY=FALSE,
  USE.NAMES=TRUE,
  mc.cores=cores
)

# make temp output dir for BOLD data
tmp.path <- paste0(
  "bold_qc_",
  paste(
    month(ymd(Sys.Date()), label=TRUE),
    year(ymd(Sys.Date())),
    sep="-"
  )
)
if(!dir.exists(here("temp", tmp.path))){
  dir.create(here("temp", tmp.path))
}

# run phylo for BOLD data
trs.list <- mapply(
  function(x,y) phylogenize(
    dir=tmp.path, 
    fas=x, 
    prefix=paste0("bold_", y), 
    verbose=opt$verbose
  ), 
  reflibs.fas, 
  prefixes, 
  SIMPLIFY=FALSE,
  USE.NAMES=TRUE
)

# plot the trees in a temp dir
mcmapply(
  function(x,y,z) plot_trees(
    tr=x, 
    df=y, 
    prefix=paste0("bold_", z)
  ), 
  trs.list, 
  reflibs.haps, 
  prefixes, 
  SIMPLIFY=FALSE,
  USE.NAMES=TRUE,
  mc.cores=cores
)

end_time <- Sys.time()
end_time-start_time

writeLines("\nPhylogenetic trees have been generated for BOLD data, PDFs are available in 'reports/bold_qc_MONTH-YEAR'")