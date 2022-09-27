
Clean_Names <- function(names, show.progress = FALSE, return_gen_sps = TRUE){
    res <- c()
    for(x in 1:length(names)){
        if(show.progress){
            cat("\r",x,"from",length(names))
        }
        if(is.na(names[x])){
            # do nothing
        }
        else{
            # replace artifacts from encoding issue
            tmp. <- gsub('<U\\+FB01>', 'fi', names[x])
            tmp. <- gsub('ÿ', '', tmp.)
            # dot within name
            tmp. <- strsplit(tmp.,"[.]")[[1]]
            tmp. = paste(tmp., collapse = " ")
            # reduces repeated white space
            tmp. <- str_squish(tmp.)
            # remove letters with accent
            tmp. <- iconv(tmp., "latin1", "ASCII", "")
            # remove parentheses
            tmp. <- strsplit(tmp., " ")[[1]]
            if(any(grepl("[(]", tmp.))){
                tmp. = tmp.[-which(grepl("[(]", tmp.))]
            }
            # remove bar
            if(any(grepl("[/]", tmp.))){
                pos <- grep("[/]", tmp.)
                tmp.[pos] = gsub("[/]"," ",tmp.[pos])
                tmp. = paste(tmp., collapse = " ")
                tmp. <- strsplit(tmp., " ")[[1]]
            }
            # Genus uppercase and species lowercase
            tmp.1 = str_to_title(tmp.[1])
            if(length(tmp.)>1){
                tmp.2 <- sapply(tmp.[-1], function(x) str_to_lower(x))
                tmp.2 = paste(tmp.1, paste(tmp.2, collapse = " "), collapse = " ")
                tmp. <- strsplit(tmp.2, " ")[[1]]
            } else {
                tmp. <- tmp.1
            }
            # fix sp. var. subspecies
            tmp.[which(tmp. == "sp")] <- "sp."
            tmp.[which(tmp. == "spp")] <- "sp."
            tmp.[which(tmp. == "spp.")] <- "sp."
            tmp.[which(tmp. == "var")] <- "var."
            tmp.[which(tmp. == "ssp")] <- "subsp."
            tmp.[which(tmp. == "ssp.")] <- "subsp."
            tmp.[which(tmp. == "subsp")] <- "subsp."
            tmp. <- gsub("[..]",".",tmp.)
            if(return_gen_sps) {
                # remove var subsp cf
                if(any(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")){
                    tmp. = tmp.[-which(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")]
                }
                if(length(tmp.)>2){
                    tmp. <- tmp.[1:2]
                }
            }
            if(length(tmp.)==1){
                tmp. = paste(tmp., "sp.", collapse = " ")
            } else {
                tmp. = paste(tmp., collapse = " ")
            }
        }
        
        res = c(res,tmp.)
    }
    return(res)
}


####################
# spnames = togo$species
# db = "ncbi"
# suggest_names = FALSE
# return_accepted_only = TRUE
# # 
# test = Find_Names(togo$species, db = "gbif", suggest_names = FALSE, return_accepted_only = TRUE)

# get_n_occ = whereas to retrieve N occurrences from the taxa. Only works if db = 'gbif'.
# basisOfRecord = for extracting gbif data. See rgbif::occ_search. Only works if get_n_occ is TRUE.
# year = temporal range of gbif occurrence records. See rgbif::occ_search. Only works if get_n_occ is TRUE.

Find_Names <- function(spnames, 
                       db = "gbif", 
                       suggest_names = FALSE, 
                       return_accepted_only = FALSE,
                       parallel = F,
                       ncores = detectCores()-2){
    cat("Searching names for",length(spnames),"taxa\n")
    
    sps_names<-query_names(sci_name = spnames, 
                           db = db,
                           suggest_names = suggest_names,
                           parallel = parallel,
                           ncores = ncores) 
    sps_names$scientificName <- paste(sps_names$genus,sps_names$specificEpithet)
    sps_names$scientificName[which(sps_names$scientificName=="NA NA")] <- NA
    if(any(duplicated(sps_names$original_search))){
        # keep the one with info
        dups <- unique(sps_names$original_search[which(duplicated(sps_names$original_search))])
        dups <- which(sps_names$original_search %in% dups)
        tokeep <- sps_names[dups,]
        tokeep <- tokeep %>%
            arrange(original_search, is.na(taxonID)) %>% 
            distinct(original_search, .keep_all = TRUE)
        sps_names <- sps_names[-dups,]
        sps_names <- rbind(sps_names,tokeep)
    }
    # all(spnames %in% sps_names$original_search)
    
    found <- length(unique(sps_names$original_search[
        grep("accepted",sps_names$taxonomicStatus)]))
    cat("Found accepted names for",found,"out of",length(spnames),"taxa\n")
    
    ## find accepted names for synonyms from accepted id
    cat("Searching accepted names for synonims from accepted id\n")
    ## Method 1 - from Accepted ID
    pos <- which(!sps_names$taxonomicStatus == "accepted" & !is.na(sps_names$acceptedNameUsageID))
    if(any(pos)){
        cat("Method 1 - Get names from accepted ID\n")
        cat("Searching for",length(pos),"taxa\n")
        togosyno <- sps_names[pos,]
        ## find names from ID
        togosyno_accepted <- find_names_from_id(id = togosyno$acceptedNameUsageID,
                                                name = togosyno$original_search,
                                                db = db)
        # feed
        togosyno_accepted <- togosyno_accepted[,which(colnames(togosyno_accepted) %in% colnames(sps_names))]
        rem <- which(sps_names$original_search %in% togosyno_accepted$original_search)
        if(any(rem)){
            sps_names <- sps_names[-rem,]
        }
        sps_names <- bind_rows(sps_names,
                               togosyno_accepted)
        
        found <- length(unique(togosyno_accepted$original_search[
            grep("accepted",togosyno_accepted$taxonomicStatus)]))
        
        cat("Found accepted names for",found,"taxa\n")
    }
    ## Method 2 - using function taxadb::filter_name
    pos <- which(!sps_names$notes == "accepted")
    if(any(pos)){
        cat("Method 2 - Filter names \n")
        cat("Searching for",length(pos),"taxa\n")
        togosyno <- sps_names[pos,]
        togosyno_accepted <- find_names_filter(x = togosyno$original_search,
                                               db = db)
        # feed
        togosyno_accepted <- togosyno_accepted[,which(colnames(togosyno_accepted) %in% colnames(sps_names))]
        rem <- which(sps_names$original_search %in% togosyno_accepted$original_search)
        if(any(rem)){
            sps_names <- sps_names[-rem,]
        }
        sps_names <- bind_rows(sps_names,
                               togosyno_accepted)
        
        found <- length(unique(togosyno_accepted$original_search[
            grep("accepted",togosyno_accepted$taxonomicStatus)]))
        
        cat("Found accepted names for",found,"taxa\n")
    }
    
    cat("--- Summary ---\n",
        "N taxa:",nrow(sps_names),"\n",
        "N taxa found:",length(grep("accepted",sps_names$taxonomicStatus)), "\n",
        "N taxa not found:", nrow(sps_names)-length(grep("accepted",sps_names$taxonomicStatus)))
    if(return_accepted_only){
        # select only accepted names
        sps_names <- sps_names[grep("accepted",sps_names$taxonomicStatus),]
    } 
    return(sps_names)
}

####################
find_names_from_id <- function(id, name, db = "gbif"){
    togosyno_accepted <- taxadb::filter_id(id, 
                                           provider = db, 
                                           type = "acceptedNameUsageID")
    name = data.frame(original_search = name, input = id)
    togosyno_accepted <- merge(togosyno_accepted, name, by = "input")
    
    if(any(duplicated(togosyno_accepted$acceptedNameUsageID))){
        # use the most similar name
        dups <- unique(togosyno_accepted$acceptedNameUsageID[which(duplicated(togosyno_accepted$acceptedNameUsageID))])
        found <- lapply(dups, function(x){
            tmp <- togosyno_accepted[which(togosyno_accepted$acceptedNameUsageID == x),]
            tmp[order(adist(unique(tmp$original_search),tmp$scientificName))[1],]
        })
        found <- do.call(rbind,found)
        togosyno_accepted <- togosyno_accepted[-which(togosyno_accepted$acceptedNameUsageID %in% dups),]
        togosyno_accepted <- rbind(togosyno_accepted,found)
    }
    return(togosyno_accepted)
}

####################
query_names <- function(sci_name,
                        replace_synonyms = TRUE,
                        suggest_names = TRUE,
                        suggestion_distance = 0.9,
                        db = "gbif",
                        rank_name = NULL,
                        rank = NULL,
                        parallel = FALSE,
                        ncores = 2,
                        export_accepted = FALSE){
    
    ori_names <- sci_name # save original names
    if(db == 'gbif'){ # gbif can find anything with a '-'
        sci_name = gsub("-"," ",sci_name)
    }
    ori_names <- data.frame(original_search = ori_names, sci_name)
    
    sps_names<-bdc_query_names_taxadb(sci_name = ori_names$sci_name,
                                      replace_synonyms = replace_synonyms,
                                      suggest_names = suggest_names,
                                      suggestion_distance = suggestion_distance,
                                      db = db,
                                      parallel = parallel,
                                      ncores = ncores,
                                      export_accepted = export_accepted) 
    sps_names$scientificName <- paste(sps_names$genus,sps_names$specificEpithet)
    sps_names$scientificName[which(sps_names$scientificName=="NA NA")] <- NA
    if(db == 'gbif'){ # return original names to data
        sps_names$id <- sps_names$original_search
        sps_names <- sps_names[,-which(colnames(sps_names)=="original_search")]
        sps_names <- merge(sps_names, ori_names, by.x = "id",by.y = 'sci_name',all.x = T)
        sps_names <- data.frame(sps_names)
        sps_names <- sps_names[,-which(colnames(sps_names)=="id")]
    } 
    return(sps_names)
}

####################
# if(db == 'gbif'){
#     # select the one with the greatest N occurrence
#     ids <- gsub("GBIF:","",togosyno_accepted$acceptedNameUsageID)
#     nocc <- get_n_occur_gbif(gbif_code = ids)
#     dups_names = data.frame(original_search=togosyno_accepted$original_search,
#                             nocc)
#     tolook <- unique(dups_names$original_search)
#     found <- lapply(tolook, function(x){
#         tmp <- dups_names[which(dups_names == x),]
#         tmp[which(tmp$nocc == max(tmp$nocc)),]
#     })
#     togosyno_accepted <- do.call(rbind,found)
# }

####################
# x=togosyno$original_search
# db="itis"

# "Alchemilla arvensis"
# "Perezia microcephala"

find_names_filter <- function(x, db = db){
    togosyno_accepted <- taxadb::filter_name(x, 
                                             provider = db)
    togosyno_accepted$original_search <- togosyno_accepted$input
    if(any(is.na(togosyno_accepted$acceptedNameUsageID))){
        togosyno_accepted <- togosyno_accepted[-which(is.na(togosyno_accepted$acceptedNameUsageID)),]
    }
    return(togosyno_accepted)
}

