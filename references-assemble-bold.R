#!/usr/bin/env Rscript

# R script to make reference databases using BOLD data and HMM processing
# Modified to use BOLD data only

## Load functions and libs
source(here::here("scripts/load-libs.R"))

# load synonyms from rfishbase
source(here("scripts/load-synonyms.R"))

# get args
option_list <- list(
  make_option(c("-t","--threads"), type="numeric"),
  make_option(c("-m","--metabarcode"), type="character")
)

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
# if running line-by-line
#opt <- NULL
#opt$threads <- 1
#opt$metabarcode <- "coi.lerayxt"

# set cores
cores <- opt$threads

## Data
# load up the species table
species.table <- read_csv(file=here("assets/species-table.csv"),show_col_types=FALSE)
# load the BOLD dump
bold.red <- read_csv(file=here("temp/bold-dump.csv"), guess_max=100000,show_col_types=FALSE)
# load up stats
stats <- read_csv(file=here("reports/stats.csv"),show_col_types=FALSE)

## Extract the frag of interest using the HMMs in hmmer
# get list of metabarodes
prefixes.list <- c("coi.lerayxt","coi.ward","12s.miya","12s.riaz","12s.valentini","12s.taberlet","16s.berry","cytb.minamoto","16s.kitano")
# split the input
prefixes.chosen <- unlist(str_split(opt$metabarcode,","))

# choose metabarcode
if(opt$metabarcode == "all") {
  prefixes.all <- c("coi.lerayxt.noprimers","coi.ward.noprimers","12s.miya.noprimers","12s.riaz.noprimers","12s.valentini.noprimers","12s.taberlet.noprimers","16s.berry.noprimers","cytb.minamoto.noprimers","16s.kitano.noprimers")
} else if (all(prefixes.chosen %in% prefixes.list)) {
  prefixes.all <- paste(prefixes.chosen,"noprimers",sep=".")
} else stop(writeLines("'-m' value must be metabarcode(s) listed in Table 1, and separated by a comma, e.g. '12s.miya,coi.ward'."))

# run hmmer
writeLines("\nExtracting metabarcode fragments with HMMER (may take several minutes) ...")
# use single thread because easier on the RAM
dat.frag.all <- lapply(prefixes.all, function(x) run_hmmer3(dir="temp", infile="mtdna-dump.fas", prefix=x, evalue="10", coords="env"))
writeLines("\nDone")

# concatenate all
dat.frag.cat <- do.call(c,dat.frag.all)

# get unique names
dat.frag.names <- unique(labels(dat.frag.cat))

# get BOLD sequences
bold_sequences <- dat.frag.names[dat.frag.names %in% bold.red$processid]

# clean BOLD data
dbs.merged.all <- bold.red %>% 
  filter(processid %in% bold_sequences) %>%
  filter(!is.na(species_name)) %>%
  mutate(source="BOLD",
         nucleotides=str_to_lower(nucleotides), 
         length=as.character(str_length(nucleotides))) %>% 
  select(source, processid, species_name, lat, lon, country, 
         institution_storing, catalognum, nucleotides, length) %>%
  rename(dbid=processid,sciNameOrig=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)


# name each DNAbin object
names(dat.frag.all) <- prefixes.all

# extract nucleotides from DNAbin objects
dat.frag.flat <- lapply(dat.frag.all, function(x) 
  mcmapply(str_flatten, as.character(x), mc.cores=cores, SIMPLIFY=TRUE, USE.NAMES=TRUE))

# convert to dataframes
dat.frag.df <- lapply(dat.frag.flat, function(x) 
  tibble(names=names(x), seqs=unlist(x), lengthFrag=str_length(seqs)))

# rename dataframes with fragment names
dat.frag.df <- mapply(function(x,y,z) 
  dplyr::rename(x, dbid=names, !!y:=seqs, !!z:=lengthFrag), 
  dat.frag.df, 
  paste("nucleotidesFrag", names(dat.frag.df), sep="."), 
  paste("lengthFrag", names(dat.frag.df), sep="."), 
  SIMPLIFY=FALSE)

# merge all dataframes
dat.frag.merged <- dat.frag.df %>% 
  purrr::reduce(full_join, by="dbid") 
 # rename(matchCol=dbid)
# join with metadata
dbs.merged.all <- dplyr::left_join(dbs.merged.all, dat.frag.merged, by="dbid")

## Add fishbase taxonomy
# clean species names
dbs.merged.all %<>% 
  mutate(sciNameBinomen=sciNameOrig,
         sciNameBinomen=str_replace_all(sciNameBinomen," sp\\. "," sp."),
         sciNameBinomen=str_replace_all(sciNameBinomen," cf\\. "," cf."),
         sciNameBinomen=str_replace_all(sciNameBinomen," aff\\. "," aff.")) %>% 
  mutate(sciNameBinomen=apply(str_split_fixed(sciNameBinomen, " ", 3)[,1:2], 1, paste, collapse=" "))

# prepare valid species reference
uk.species.valid <- species.table %>% 
  distinct(fbSpecCode, validName, class, order, family, genus, commonName) %>% 
  mutate(rank=if_else(grepl(" ",validName),"species","genus"))
uk.species.genera <- uk.species.valid %>% 
  filter(rank=="genus") %>% 
  pull(validName)

# annotate with fishbase data
dbs.merged.all %<>% 
  mutate(fbSpecCode=pull(fishbase.synonyms.acc,SpecCode)[match(sciNameBinomen,pull(fishbase.synonyms.acc,synonym))]) %>% 
  mutate(fbSpecCode=if_else(is.na(fbSpecCode),
                            pull(fishbase.synonyms.syn,SpecCode)[match(sciNameBinomen,pull(fishbase.synonyms.syn,synonym))],
                            fbSpecCode)) %>%
  mutate(genus=str_split_fixed(sciNameBinomen," ",2)[,1]) %>%
  mutate(rank=if_else(genus %in% uk.species.genera,"genus","species")) %>%
  mutate(sciNameValid=if_else(rank=="species",
                              pull(uk.species.valid,validName)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))],
                              sciNameBinomen)) %>%
  mutate(genus=if_else(rank=="species",str_split_fixed(sciNameValid," ",2)[,1],genus))

# handle missing taxa
missing <- dbs.merged.all %>% 
  filter(is.na(sciNameValid)) %>% 
  pull(sciNameOrig) %>% 
  unique()
if(length(missing)>0) {
  writeLines(paste("\nThe following taxa could not be found in the species database and have been dropped:",
                   paste(missing,collapse=", ")))
  dbs.merged.all %<>% filter(!is.na(sciNameValid))
}

# report genus level taxa
only.genera <- dbs.merged.all %>% 
  filter(rank=="genus") %>% 
  distinct(sciNameValid) %>% 
  pull()
if(length(only.genera)>0) {
  writeLines(paste("\nThe following taxa were searched for at the genus level:",
                   paste(only.genera,collapse=", ")))
}

# report updated names
updated <- dbs.merged.all %>% 
  filter(sciNameOrig != sciNameValid) %>% 
  select(sciNameOrig,sciNameValid) %>% 
  arrange(sciNameOrig) %>% 
  distinct()
if(nrow(updated)>0) {
  writeLines("\nThe following taxa had their BOLD names updated using FishBase:")
  print(updated,n=Inf)
}

# add taxonomy
dbs.merged.all %<>% 
  mutate(phylum="Chordata") %>%
  left_join(distinct(uk.species.valid,class,order,family,genus),by="genus")

# check taxonomy completeness
no.tax <- dbs.merged.all %>% 
  filter(is.na(class) | is.na(order) | is.na(family)) %>% 
  select(class,order,family,sciNameOrig,sciNameValid) %>% 
  arrange(sciNameOrig) %>% 
  distinct()
if(nrow(no.tax)>0) {
  writeLines("\nThe following taxa could not be assigned taxonomy. Consider adding these to the species table if you want to keep them.")
  print(no.tax,n=Inf)
}

# prepare final output
dbs.merged.info <- dbs.merged.all %>% 
  select(-matches("Frag")) %>% 
  select(source, dbid, sciNameValid, phylum, class, order, family, genus, sciNameOrig, fbSpecCode,
         country, catalogNumber, institutionCode, decimalLatitude, decimalLongitude, length, nucleotides)

dbs.merged.seqs <- dbs.merged.all %>% 
  select(matches("Frag|dbid"))

dbs.merged.final <- left_join(dbs.merged.info, dbs.merged.seqs, by="dbid") %>%
  arrange(class, order, family, genus, sciNameValid) %>% 
  filter(!is.na(nucleotides))

# write output
# write out a gzipped file (orig is too big for github)
writeLines("\nWriting out reference library to 'assets/reference-library-master.csv.gz' ...")
write_csv(dbs.merged.final, file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
writeLines("\nAll operations completed!\nPlease read previous messages in case of error")