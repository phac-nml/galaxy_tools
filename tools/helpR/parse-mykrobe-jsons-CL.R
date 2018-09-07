# Parsing JSONs from mykrobe to replace Dan's extreme workflow
# Take the json output from mykrobe, rearrange, output for LIMS
# Adrian Zetner
# August 2018

# Libraries ####
library(jsonlite)
library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(optparse)

# Define custom functions, variables, and paths. Collect and use CL arguments ####

# Here's a function to recreate that output table from the input JSON files

getResults <- function(listelement){
  # Define list levels for various things
  species <- names(listelement[[1]][["phylogenetics"]][["species"]]) 
  lineage <- names(listelement[[1]][["phylogenetics"]][["lineage"]])
  phylo_group <- names(listelement[[1]][["phylogenetics"]][["phylo_group"]])
  if("Non_tuberculosis_mycobacterium_complex" %in% phylo_group){
    warning(paste("Non-tuberculosis mycobacteria detected in file ", names(listelement), ". Skipping.", sep = ""))
    return()}
    
  # Start building a list of all your various elements
  temp <- list(mykrobe_version = listelement[[1]][["version"]][["mykrobe-predictor"]],
               file = names(listelement), # One element
               plate_name = "test", # This probably needs changing
               sample = "sequence_calls", # Likewise change this
               phylo_group = phylo_group, # As above
               species = species, # As above
               lineage = lineage, # As above
               # The following expressions drill down into the list elements and pull out what is needed. 
               # It's inelegant and vulnerable to changes in the input formats but if they're consistent it'll work
               phylo_group_per_covg = listelement[[1]][["phylogenetics"]][["phylo_group"]][[phylo_group]][["percent_coverage"]],
               species_per_covg = listelement[[1]][["phylogenetics"]][["species"]][[species]][["percent_coverage"]],
               lineage_per_covg = listelement[[1]][["phylogenetics"]][["lineage"]][[lineage]][["percent_coverage"]],
               phylo_group_depth = listelement[[1]][["phylogenetics"]][["phylo_group"]][[phylo_group]][["median_depth"]],
               species_depth = listelement[[1]][["phylogenetics"]][["species"]][[species]][["median_depth"]],
               lineage_depth = listelement[[1]][["phylogenetics"]][["lineage"]][[lineage]][["median_depth"]],
               Mykrobe_Resistance_probe_set = basename(listelement[[1]][["probe_sets"]][2]) # Is it always the second?
  )
  
  # Super cool nested and vectorized (for SPEED!) functions to grab the predictions for drug sensitivity and gene variants
  # Both produce character vectors of the same length as the number of drugs tested in the same order
  # All of these also check if there are missing values in drug/susceptibility/variant elements and adds the column anyhow
  
  if(length(map_chr(listelement[[1]][["susceptibility"]],  "predict")) != 0){
    temp$susceptibility <- map_chr(listelement[[1]][["susceptibility"]],  "predict")  
  }else{
    temp$susceptibility <- NA
  }
  
  if(length(names(listelement[[1]][["susceptibility"]])) != 0){
    temp$drug <- names(listelement[[1]][["susceptibility"]])  
  }else{
    temp$drug <- NA
  }
  
  mapped.variants <- map(listelement[[1]][["susceptibility"]], # Dig into the lists, pull out variants and collapse into chr vector
                       ~ imap(.x[["called_by"]],  # imap is shorthand for map2(x, names(x), ...), calling .y gets you the name / index of the current element
                              ~ paste(.y, 
                                      .x[["info"]][["coverage"]][["alternate"]][["median_depth"]],
                                      .x[["info"]][["coverage"]][["reference"]][["median_depth"]],
                                      .x[["info"]][["conf"]],
                                      sep = ":"))) %>% 
    map_chr(~ paste(.x, collapse = "__"))
  
  if(length(mapped.variants) != 0){
    temp$`variants (gene:alt_depth:wt_depth:conf)` <- mapped.variants
  }else{
    temp$`variants (gene:alt_depth:wt_depth:conf)` <- NA
  }
  
  temp$`genes (prot_mut-ref_mut:percent_covg:depth)` <- NA # This is dumb but it makes it match that last column. What is it there for?
  
  # Take that list and mash all the elements together as columns in a tibble, recycling as needed to fill in space
  # eg. phylo_group is repeated/recycled as many times as there are drugs tested
  as_tibble(temp)
}

# Get command line arguments with optparse
option_list = list(
  make_option(c("-f", "--file"), 
              type="character", 
              default=NULL, 
              help="dataset file name or comma separated names: eg. file1,file2,file3", 
              metavar="character"),
  make_option(c("-d", "--dir"), 
              type="character", 
              default=NULL, 
              help="directory of json files", 
              metavar="character"),
  make_option(c("-v", "--version"), 
              type="character", 
              default="", 
              help="Mykrobe Workflow Version", 
              metavar="character"),
  make_option(c("-D", "--depth"), 
              type="integer", 
              default=5, 
              help="Minimum depth of coverage [default= %default]", 
              metavar="integer"),
  make_option(c("-c", "--conf"), 
              type="integer", 
              default=10, 
              help="Minimum genotype confidence for variant genotyping [default= %default]", 
              metavar="integer"),
  make_option(c("-n", "--name"), 
              type="character", 
              default="", 
              help="Name of the run", 
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file) & is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied to input file or directory", call.=FALSE)
}

# Parameters to take from Galaxy as args or however works best
params <- c("",  # Lims_Comment 
            "",  # Lims_INTComment 
            opt$version,  # Mykrobe_Workflow_Version
            opt$depth,  # Mykrobe_min_depth_default_5
            opt$conf,  # Mykrobe_min_conf_default_10
            "",                             # LIMS_file - empty as it's an upload field in LIMS
            opt$name)  # LIMS_filename

names(params) <- c("Lims_Comment", 
                   "Lims_INTComment",
                   "Mykrobe_Workflow_Version",
                   "Mykrobe_min_depth_default_5",
                   "Mykrobe_min_conf_default_10", 
                   "LIMS_file", 
                   "LIMS_filename")


# A default report in the order LIMS requires

# Make a default dataframe to combine the rest into and enforce column order / fill missing ones with NAs
columns <- c("file",
             "Mykrobe_fabG1",
             "Mykrobe_katG",
             "Mykrobe_ahpC",
             "Mykrobe_inhA",
             "Mykrobe_ndh",
             "Isoniazid_R_mutations",
             "Isoniazid_Prediction",
             "Mykrobe_rpoB",
             "Rifampicin_R_mutations",
             "Rifampicin_Prediction",
             "Mykrobe_embB",
             "Mykrobe_embA",
             "Ethambutol_R_mutations",
             "Ethambutol_Prediction",
             "Mykrobe_pncA",
             "Mykrobe_rpsA",
             "Pyrazinamide_R_mutations",
             "Pyrazinamide_Prediction",
             "Mykrobe_gyrA",
             "Moxifloxacin_Ofloxacin_R_mutations",
             "Moxifloxacin_Ofloxacin_Prediction",
             "Mykrobe_rpsL",
             "Mykrobe_Streptomycin_rrs",
             "Mykrobe_Streptomycin_gid",
             "Streptomycin_R_mutations",
             "Streptomycin_Prediction",
             "Mykrobe_Amikacin_rrs",
             "Amikacin_R_mutations",
             "Amikacin_Prediction",
             "Mykrobe_Capreomycin_rrs",
             "Mykrobe_Capreomycin_tlyA",
             "Capreomycin_R_mutations",
             "Capreomycin_Prediction",
             "Mykrobe_Kanamycin_rrs",
             "Mykrobe_Kanamycin_eis",
             "Kanamycin_R_mutations",
             "Kanamycin_Prediction",
             "Lims_Comment",
             "Lims_INTComment",
             "Mykrobe_Workflow_Version",
             "mykrobe_version",
             "Mykrobe_Resistance_probe_set",
             "Mykrobe_min_depth_default_5",
             "Mykrobe_min_conf_default_10",
             "LIMS_file",
             "LIMS_filename")

report <- setNames(data.frame(matrix("", ncol = length(columns), nrow = 1), stringsAsFactors = F), columns)


# List of drugs that are tested
all_drugs <- c("Isoniazid", 
               "Rifampicin", 
               "Ethambutol", 
               "Pyrazinamide", 
               "Moxifloxacin_Ofloxacin", 
               "Streptomycin",
               "Amikacin",
               "Capreomycin",
               "Kanamycin")

# Do Stuff ####

# Import all the JSON files into a list of lists format ####

# A directory of files

# in.dir <- c("data/mykrobe_predictor_files_for_adrian/Data - Json files from Mykrobe/", "data/failures/")
if (is.null(opt$file)){
  # opt$dir is used to get the list of files, a vector of non-duplicated files is then passed to map
  files <- list.files(path = opt$dir, 
                      pattern = "*.json", # HERE PHIL! CHANGE THIS TO WHAT YOU NEED IT TO BE IN GALAXY
                      full.names = T)
}else{
  files <- unlist(strsplit(opt$file, ","))
}

files <- files[!duplicated(basename(files))]

list.of.json.files <- map(files, 
                          ~ fromJSON(.x, simplifyDataFrame = F)
)


# Apply that getResults function to each element in your list then bash it together into a final report

temp <- map(list.of.json.files, getResults) %>% 
  bind_rows()


# Predictions of resistance or susceptibility
predictions.table <- 
  temp %>%
  select(file, drug, susceptibility) %>% 
  mutate(drug = paste(drug, "_Prediction", sep = "")) %>% 
  spread(drug, susceptibility, fill = "failed") %>% 
  select(-starts_with("NA"))

# Variants
# Multiple resistance mutations and confidence per drug in the X_R_mutations column
# Actual protein changes in Mykrobe_X columns

variants.temp <- 
  temp %>% 
  select(file, drug, variants = `variants (gene:alt_depth:wt_depth:conf)`) %>% 
  # bind_rows(tibble(drug = all_drugs)) %>% # This was to add the rest of the drugs, but why bother if just populating a table later?
  # complete(file, drug) %>% 
  # filter(!is.na(file)) %>%
  # filter(!is.na(drug)) %>%
  mutate(variants = replace(variants, variants == "", NA)) %>% # Make missing data consistent
  filter(!is.na(variants)) %>% # Then get rid of it
  mutate(tempcols = paste(drug, "R_mutations", sep = "_")) %>% 
  mutate(R_mutations = variants) %>% 
  mutate(variants = strsplit(variants, "__")) %>% # Split the mutations across rows (list first then split across rows)
  unnest(variants) %>% 
  separate(variants, c("gene", "mutation"), "_") %>% 
  mutate(columnname = ifelse(gene %in% c("tlyA", "rrs", "gid"), # Check for columns that include the drug name or not and paste accordingly
                             paste("Mykrobe", drug, gene, sep = "_"),
                             paste("Mykrobe", gene, sep = "_"))) %>% 
  # Extract out the mutation information with a regex that covers all potential genes
  # This regex looks for whatever is ahead of the first colon and after the last hyphen
  mutate(mutation = str_match(mutation, "(.*)-.*:")[,2]) %>%
  select(file, tempcols, R_mutations, columnname, mutation)

# Split each kind of variants into its own temp table then merge
variants.1 <- 
  variants.temp %>% 
  select(file, tempcols, R_mutations) %>% 
  distinct() %>% 
  spread(tempcols, R_mutations)

variants.2 <- 
  variants.temp %>% 
  select(file, columnname, mutation) %>% 
  spread(columnname, mutation)

variants.table <- full_join(variants.1, variants.2, by = "file")

# Make a report ####

report <- 
  temp %>% 
  select(file, mykrobe_version, Mykrobe_Resistance_probe_set) %>% # Get important info from initial table
  distinct() %>% # Drop duped rows and combine all the tables together
  full_join(variants.table) %>% 
  full_join(predictions.table) %>% 
  bind_rows(report) %>% # Use bind_rows to add columns (eg. unteseted drugs) to the final output 
  filter(file != "")

# Only add the 'no mutation' replacement to the columns that actually have a result
report <- 
  report %>%
  filter_at(vars(ends_with("_Prediction")), any_vars(. != "failed")) %>% 
  mutate_at(vars(starts_with("Mykrobe_")), funs(replace(., is.na(.), "No Mutation"))) %>% 
  full_join(anti_join(report, ., by = "file")) %>% 
  select(columns)

# Add in the parameters fed from Galaxy using named character vector
report <- 
  report %>% 
  mutate(
    Lims_Comment = params["Lims_Comment"],
    Lims_INTComment = params["Lims_INTComment"],
    Mykrobe_Workflow_Version = params["Mykrobe_Workflow_Version"],
    Mykrobe_min_depth_default_5 = params["Mykrobe_min_depth_default_5"],
    Mykrobe_min_conf_default_10 = params["Mykrobe_min_conf_default_10"],
    LIMS_file = params["LIMS_file"],
    LIMS_filename = params["LIMS_filename"]
  )

#View(report)

# Write some output
# Report as is
write.csv(report, "data/output-report.txt", row.names = F)

# Select specific columns from temp and output them
temp %>% 
  select(file, 
         phylo_group, 
         species, 
         lineage, 
         phylo_group_per_covg, 
         species_per_covg, 
         lineage_per_covg, 
         phylo_group_depth, 
         species_depth, 
         lineage_depth) %>%
  distinct() %>%
  write.csv("data/output-jsondata.txt", row.names = F)
