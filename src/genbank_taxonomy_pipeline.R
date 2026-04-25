
#an R script to cross reference each GenBank accession ID hit and 
#retrieve the Organism level information from NCBI


####################################################################
####################################################################
#1. LOADING THE REQUIRED PACKAGES
#and installing if missing any needful

required_packages <- c("tidyverse", "rentrez", "readxl")

install_if_missing <- function(pkg)  {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#applying the function to all required packages 
invisible(lapply(required_packages, install_if_missing))


###################################################################
#2. SETTING WORKING DIRECTORY

getwd()
setwd("/Users/nickm/msc_bidpg_2026/bioinformatics_module/biological_databases_ontology")
getwd()


####################################################################
####################################################################
#3. LOADING AND PREPING INPUT FILES
#with error handling

#define input file, and
#path containing the input file
filetype <- "xlsx"

filepath <- "/Users/nickm/msc_bidpg_2026/bioinformatics_module/biological_databases_ontology"

#list files in the directory, 
#return full paths
if (!dir.exists(filepath)) {
  stop("Error: Specified directory does not exist.")
}

if (grepl("csv", filetype, ignore.case = TRUE)) {
  files <- list.files(filepath, pattern = "\\.csv$", full.names = TRUE)
} else if (grepl("tsv", filetype, ignore.case = TRUE)) {
  files <- list.files(filepath, pattern = "\\.txt$", full.names = TRUE)
} else if (grepl("xlsx", filetype, ignore.case = TRUE)) {
  files <- list.files(filepath, pattern = "\\.xlsx$", full.names = TRUE)
} else {
  stop("Error: Unsupported file type. Use csv, tsv, or xlsx.")
}

#stop execution if there are no files found
if (length(files) == 0) {
  stop("Error: No input files found in the specified directory.")
}

#set wd to file path, and
#track total number of processed rows, then
#process the start time
setwd(filepath)
nrows_parsed <- 0
start <- Sys.time()


#################################################################################
##################################################################################
#4. BUILDING THE NCBI QUERY FUNCTION

#store already queried IDs,
accession_cache <- list()

#function to query NCBI, and
#extract organism info
getOrganismInfo <- function(accession_id) {
  
  record <- tryCatch({
    rentrez::entrez_fetch(db = "protein", id = accession_id, rettype = "text")
  }, error = function(e_protein) {
    tryCatch({
      rentrez::entrez_fetch(db = "nuccore", id = accession_id, rettype = "text")
    }, error = function(e_nuccore) {
      return(list(ID_host = NA, ID_country = NA, ID_organism = NA))
    })
  })
  
  #split text into individual lines
  info_lines <- strsplit(record, "\n")[[1]]
  
  #then extract host information,
  #and country info if available
  host_line <- grep("/host=", info_lines, ignore.case = TRUE, value = TRUE)
  country_line <- grep("/country=", info_lines, ignore.case = TRUE, value = TRUE)
  
  #if value exists, record the first match
  host_value <- ifelse(length(host_line) > 0, host_line[1], NA)
  country_value <- ifelse(length(country_line) > 0, country_line[1], NA)
  
  # ORGANISM extraction to get precise match from record
  org_idx <- grep("^\\s*ORGANISM", info_lines)
  ref_idx <- grep("REFERENCE", info_lines)
  
  if (length(org_idx) == 0 || length(ref_idx) == 0) {
    organism_value <- NA
  } else {
    organism_value <- paste(info_lines[org_idx:(ref_idx[1] - 1)], collapse = "")
  }
  
  #return results as named list
  return(list(ID_host = host_value, ID_country = country_value, ID_organism = organism_value))
}

#wrapper function,
#if accession already queried, return stored result
#otherwise query NCBI, and
#store results in cache
getOrganismInfoCached <- function(accession_id) {
  
  if (accession_id %in% names(accession_cache)) {
    return(accession_cache[[accession_id]])
  }
  
  result <- getOrganismInfo(accession_id)
  accession_cache[[accession_id]] <<- result
  
  Sys.sleep(0.4)
  
  return(result)
}


##################################################################################
#################################################################################
#5. MAIN PIPELINE LOOP
#reads, enriches with more information, then saves the output in the directory

for (i in files) {
  
  file_name <- basename(i)
  
  input <- tryCatch({
    if (grepl("csv", filetype, ignore.case = TRUE)) {
      read.delim(i, sep = ",", header = TRUE)
    } else if (grepl("tsv", filetype, ignore.case = TRUE)) {
      read.delim(i, sep = "\t", header = TRUE)
    } else {
      readxl::read_excel(i)
    }
  }, error = function(e) {
    message("Error reading file: ", file_name)
    return(NULL)
  })
  
  if (is.null(input)) next
  
  #ensure the required colum exists
  if (!"sseqid" %in% colnames(input)) {
    message("Skipping file: ", file_name, " - 'sseqid' column not found")
    next
  }
  
  input <- input[!(is.na(input$sseqid) | input$sseqid == ""), ]
  nrows_parsed <- nrows_parsed + nrow(input)
  
  # pre-allocation
  ID_host <- rep(NA, nrow(input))
  ID_country <- rep(NA, nrow(input))
  ID_organism <- rep(NA, nrow(input))
  
  
  # enrichment loop
  for (j in seq_len(nrow(input))) {
    
    accession_id <- input$sseqid[j]
    if (is.na(accession_id) || accession_id == "") next
    
    # skip already processed accessions
    if (!is.na(ID_organism[j])) next
    
    organism_info <- tryCatch(
      getOrganismInfoCached(accession_id),
      error = function(e) list(ID_host = NA, ID_country = NA, ID_organism = NA)
    )
    
    ID_host[j] <- organism_info$ID_host
    ID_country[j] <- organism_info$ID_country
    ID_organism[j] <- organism_info$ID_organism
    
    #return message after running every 50 accessions
    if (j %% 50 == 0) message("Processed ", j, " rows in ", file_name)
  }
  
  # assigning back once for faster processing
  input$ID_host <- ID_host
  input$ID_country <- ID_country
  input$ID_organism <- ID_organism
  
  
  # cleaning
  input$ID_host <- gsub("/host=", "", input$ID_host)
  input$ID_country <- gsub("/country=", "", input$ID_country)
  
  # "safer" cleaning only on character columns
  input <- input %>%
    mutate(across(where(is.character), ~ gsub("\"", "", .))) %>%
    mutate(across(where(is.character), stringr::str_trim))
  
  
  # taxonomy splitting to columns
  taxonomy <- data.frame(input$ID_organism)
  colnames(taxonomy) <- "raw_taxonomy"
  
  taxonomy <- tidyr::separate(
    taxonomy,
    raw_taxonomy,
    into = paste0("taxonomy_", 1:20),
    sep = ";\\s*",
    fill = "right"
  )
  
  taxonomy <- taxonomy %>%
    mutate(across(where(is.character), stringr::str_trim))
  
  taxonomy <- taxonomy[, colSums(is.na(taxonomy)) < nrow(taxonomy)]
  
  output <- cbind(input, taxonomy)
  
  # file naming
  write.csv(output, paste0(tools::file_path_sans_ext(file_name), "_parsed_output.csv"), row.names = FALSE)
  
  message("Finished: ", file_name)
}


####################################################################
#6. PERFORMANCE LOGGING

#record end time,
#create performance log txt file
#write summary stats, and
#close file connection
end <- Sys.time()

performance <- file(paste0(filepath, "/performance.txt"))

writeLines(
  c(
    paste0("Start time: ", start),
    paste0("End time: ", end),
    paste0("Lines parsed: ", nrows_parsed)
  ),
  performance
)

close(performance)