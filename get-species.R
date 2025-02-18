#!/usr/bin/env Rscript
# script to make a species and synonym list for African countries
# scripts/get-species.R -s true
# load functions and libs
source(here::here("scripts/load-libs.R"))

# load synonyms from rfishbase
source(here("scripts/load-synonyms.R"))

# load fishbase tables
fishbase.country <- suppressMessages(rfishbase::country(server="fishbase",version="latest"))
fishbase.taxonomy <- suppressMessages(rfishbase::load_taxa(server="fishbase",version="latest"))
fishbase.species <- suppressMessages(rfishbase::species(server="fishbase",version="latest"))

# get args
option_list <- list( 
  make_option(c("-s","--synonyms"), type="character", default = "true")
)

# set args
opt <- parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

# status
writeLines(paste0("\nSearching FishBase ..."))

### Species lists

# load countries
# filter on African countries - use ISO country codes
african_countries <- c("012", "024", "204", "072", "108", "120", "132", "140", "148", "174", "178", "180", "384", "262", "818", "226", "232", "231", "270", "288", "324", "624", "404", "426", "430", "450", "454", "466", "478", "480", "504", "508", "516", "562", "566", "646", "678", "686", "690", "694", "706", "710", "728", "729", "748", "834", "800", "894")

species.list <- fishbase.country %>% 
  filter(C_Code %in% african_countries) %>%
  distinct(C_Code, SpecCode, country) %>%
  left_join(distinct(fishbase.taxonomy, SpecCode, Species), by = join_by(SpecCode)) %>%
  left_join(distinct(fishbase.species, SpecCode, Fresh, Brack), by = join_by(SpecCode)) %>%
  filter(Fresh == 1)

### Synonyms

# join the countries and synonyms tables
species.list.syn <- species.list %>% 
  left_join(bind_rows(fishbase.synonyms.acc, fishbase.synonyms.syn), by = join_by(SpecCode), relationship = "many-to-many")

### Taxonomy and common names

# add the taxonomy
species.list.tax <- species.list.syn %>% 
  left_join(distinct(fishbase.taxonomy, SpecCode, Genus, Family, Order, Class), by = join_by(SpecCode))

# add the common names
species.list.com <- species.list.tax %>% 
  left_join(distinct(fishbase.species, SpecCode, FBname), by = join_by(SpecCode))

### Format

# format
species.list.form <- species.list.com %>% 
  rename(speciesName = synonym, status = Status, fbSpecCode = SpecCode, validName = Species, class = Class, order = Order, family = Family, genus = Genus, commonName = FBname) %>% 
  mutate(commonSpecies = TRUE) %>%
  relocate(speciesName, status, fbSpecCode, validName, class, order, family, genus, commonName, commonSpecies) %>% 
  arrange(class, order, family, genus, validName, status, speciesName) %>%
  select(!all_of(c("C_Code", "country", "SynCode", "TaxonLevel")))

# keep or remove synonyms
if (opt$synonyms == "true") {
  species.list.form <- species.list.form
} else if (opt$synonyms == "false") {
  species.list.form <- species.list.form %>% filter(status == "accepted name")
} else stop(writeLines("'-s' value must be 'true' or 'false'."))

### Clean up and write out

# get stats
if (opt$synonyms == "true") {
  acc <- species.list.form %>% count(status) %>% filter(status == "accepted name") %>% pull(n)
  syn <- species.list.form %>% count(status) %>% filter(status == "synonym") %>% pull(n)
} else if (opt$synonyms == "false") {
  acc <- species.list.form %>% count(status) %>% filter(status == "accepted name") %>% pull(n)
  syn <- "0"
} else stop(writeLines("'-s' value must be 'true' or 'false'."))

# write out
species.list.form %>% 
  filter(!is.na(speciesName)) %>% 
  write_csv(file = "assets/species-table.csv")

# print info
writeLines(paste0("\nWriting out species list for African countries comprising ", acc, " accepted names and ", syn, " synonyms, to 'assets/species-table.csv'.\n"))

# to give warning if NAs in speciesName
if (nrow(species.list.form) != nrow(filter(species.list.form, !is.na(speciesName)))) {
  writeLines("\nWarning! The following FishBase species names were removed as NA:")
  species.list.form %>% 
    filter(is.na(speciesName)) %>% 
    select(speciesName, fbSpecCode, validName) %>%
    print(n = Inf)
}
