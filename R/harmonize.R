## function to harmonize taxonomy 
harmonize <- function(sp_names){
  
  library(taxadb)
  library(bdc)
  source("R/clean_taxa_functions.R")
  
  sp_names = Clean_Names(sp_names,
                         return_gen_sps = TRUE)
  
  
  tofind <- data.frame(matrix(nrow = length(sp_names), ncol = 8))
  names(tofind) = c("scientificName", "kingdom", "phylum", "class", "order", "family", "db", "db_code")
  
  tofind <- data.frame(species = sp_names, tofind)
  
  tofind[,1:ncol(tofind)] = lapply(tofind[,1:ncol(tofind)], as.character) 
  
  togo <- tofind[which(is.na(tofind$scientificName)),]
  
  mycols <- c("speciesKey", "kingdom","phylum","class","order","family", "species")
  
  ## GBIF
  # retrieve sp names
  cl <- makeCluster(detectCores()-2)
  clusterExport(cl, c("togo","standardize_taxa"))
  
  gbif_names <- pblapply(togo$species, function(x){
    try(standardize_taxa(data.frame(verbatimScientificName = x), 
                         fuzzy = FALSE,
                         silent = TRUE))
  }, cl = cl)
  
  stopCluster(cl)
  
  rem <- sapply(gbif_names, class)
  rem <- which(rem=="try-error")
  
  if(length(rem) != 0) {
    gbif_names <- gbif_names[-rem]
  }
  
  gbif_names <- rbindlist(gbif_names)
  if(length(which(is.na(gbif_names$scientificName))) > 0) {
    gbif_names <- gbif_names[-which(is.na(gbif_names$scientificName)),]
  }
  gbif_names <- gbif_names[which(gbif_names$taxonRank=='species'),]
  
  # remove duplicates
  if(any(duplicated(gbif_names$verbatimScientificName))){
    gbif_names <- gbif_names[-which(duplicated(gbif_names$verbatimScientificName)),]
  }
  
  any(!gbif_names$original_search %in% sp_names)
  
  cat("--- Summary ---\n",
      "N taxa:",nrow(togo),"\n",
      "N taxa found:",nrow(gbif_names), "\n",
      "N taxa not found:", nrow(togo)-nrow(gbif_names))
  
  gbif_names <- gbif_names[,c("verbatimScientificName","scientificName","kingdom","phylum","class","order","family","taxonID")]
  names(gbif_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
  gbif_names$db <- "gbif"
  gbif_names$db_code <- gsub("http://www.gbif.org/species/","",gbif_names$db_code)
  gbif_names$db_code <- paste("GBIF:",gbif_names$db_code,sep = "")
  
  # Feed
  tofind <- tofind %>% 
    rows_patch(gbif_names, 
               by = "species")
  
  gc()
  
  ## ITIS
  togo <- tofind[which(is.na(tofind$scientificName)),]
  
  ## if no species left to find 
  if (nrow(togo) == 0) {
    return(tofind)
  }
  else {
    
    # retrieve sp names
    itis_names <- Find_Names(spnames = togo$species,
                             db = "itis",
                             suggest_names = FALSE, # using exact names
                             return_accepted_only = TRUE) 
    
    # remove duplicates
    if(any(duplicated(itis_names$original_search))){
      itis_names <- itis_names[-which(duplicated(itis_names$original_search)),]
    }
    
    any(!itis_names$original_search %in% sp_names)
    
    itis_names <- itis_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
    names(itis_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
    itis_names$db <- "itis"
    
    # Feed
    tofind <- tofind %>% 
      rows_patch(itis_names, 
                 by = "species")
    
    gc()
    
    
    ## NCBI
    togo <- tofind[which(is.na(tofind$scientificName)),]
    
    # retrieve sp names
    if(nrow(togo) == 0) {
      return(tofind)
    }
    else {
      ncbi_names <- Find_Names(spnames = togo$species,
                               db = "ncbi",
                               suggest_names = FALSE, # using exact names
                               return_accepted_only = TRUE) 
      
      # remove duplicates
      if(any(duplicated(ncbi_names$original_search))){
        ncbi_names <- ncbi_names[-which(duplicated(ncbi_names$original_search)),]
      }
      
      any(!ncbi_names$original_search %in% sp_names)
      
      ncbi_names <- ncbi_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
      names(ncbi_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
      ncbi_names$db <- "ncbi"
      
      if(nrow(ncbi_names) != 0) {
        # some names from ncbi come with duplicated genus
        # remove duplicated names
        new <- sapply(ncbi_names$scientificName, function(x){
          tmp <- strsplit(x," ")[[1]]
          if(any(duplicated(tmp))){
            paste(tmp[-1], collapse = " ")
          } else {
            x
          }
        })
        
        ncbi_names$scientificName <- new
      }
      
      # Feed
      tofind <- tofind %>% 
        rows_patch(ncbi_names, 
                   by = "species")
      
      gc()
      
      
      ## IUCN
      togo <- tofind[which(is.na(tofind$scientificName)),]
      
      if(nrow(togo) == 0) {
        return(tofind)
      }
      else {
        # retrieve sp names
        ## PS: Cant filter on iucn method
        iucn_names<-bdc_query_names_taxadb(togo$species, 
                                           db = "iucn",
                                           suggest_names = FALSE) # using exact names
        iucn_names$scientificName <- paste(iucn_names$genus,iucn_names$specificEpithet) 
        
        # remove duplicates
        if(any(duplicated(iucn_names$original_search))){
          iucn_names <- iucn_names[-which(duplicated(iucn_names$original_search)),]
        }
        
        any(!iucn_names$original_search %in% sp_names)
        
        # select only accepted names
        iucn_names <- iucn_names[which(iucn_names$taxonomicStatus == "accepted"),]
        
        cat("--- Summary ---\n",
            "N taxa:",nrow(togo),"\n",
            "N taxa found:",length(which(iucn_names$taxonomicStatus == "accepted")), "\n",
            "N taxa not found:", nrow(togo)-length(which(iucn_names$taxonomicStatus == "accepted")))
        
        iucn_names <- iucn_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
        names(iucn_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
        iucn_names$db <- "iucn"
        
        # Feed
        tofind <- tofind %>% 
          rows_patch(iucn_names, 
                     by = "species")
        
        gc()
        
        ## fix issues with Scientific names 
        new <- sapply(tofind$scientificName, function(x){
          tmp <- strsplit(x," ")[[1]]
          if(any(duplicated(tmp))){
            paste(tmp[-1], collapse = " ")
          } else {
            x
          }
        })
        
        tofind$scientificName <- new
      }
      
      ## return
      return(tofind)
    }
  }
}