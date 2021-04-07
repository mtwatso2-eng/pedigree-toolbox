isStandardFormat <- function(Code){
  grepl("^[0-9]+-[0-9]+$", Code)
}

getCrossFromStandardFormat <- function(Code){
  strsplit(Code, "-", fixed = T)[[1]][1]
}

isOrnamental <- function(germplasmName){
  grepl("^[0-9-]*$", germplasmName)
}

getParents <- function(germplasmName){
  if(isOrnamental(germplasmName))
    germplasmName <- getCrossFromStandardFormat(germplasmName)
  germplasm <- parents %>% filter(Code == germplasmName)
  if(nrow(germplasm) > 0)
    return(c(germplasm$Maternal, germplasm$Paternal))
  germplasmParents <- c(NA, NA)
  try(silent = T, {
    germplasmParents <- strsplit(brapi_get_germplasm(brapi_db()$sweetpotatobase, germplasmName = germplasmName)$pedigree, "/")[[1]]
  })
  return(germplasmParents)
}

pedigree <- local(function(germplasmName){
  recursiveGetParents <- function(germplasmName, thisPedigree = list()){
    try({thisParents <- getParents(germplasmName)})
    if(exists("thisParents")){
      for(i in 1:2){
        if(!thisParents[i] %in% c(NA, "NA", "")){
          thisPedigree[[i]] <- c(recursiveGetParents(thisParents[i]))
          names(thisPedigree)[i] <- thisParents[i]
        }
        else{
          thisPedigree[[i]] <- NA
        }
      }
    }
    return(thisPedigree)
  } 
  thisPedigree <- list()
  thisPedigree[[germplasmName]] <- recursiveGetParents(germplasmName)
  return(thisPedigree)
})

nestDepth <- function(nestedList){
  i <- 1
  while(TRUE){
    if(is.null(nestSlice(nestedList, i)))
      return(i-1)
    i <- i + 1
  }
}

nestSlice <- function(nestedList, n = 1){
  map_depth(nestedList, n - 1, function(x){names(x)}) %>% unlist %>% unname
}

getPercentParents <- function(germplasmName){
  # look for pedigree unknowns (terminal and OP parents) by...
  unknowns <- germplasmName %>%
    #...getting the pedigree of germplasmName...
    pedigree(.) %>%
    #...looking for places in the pedigree that end in an unknown (NA) parent...
    list.search(., all(is.na(.))) %>%
    #...getting the names of these unknowns...
    names(.) %>%
    #...and undo the work that list.search does to make results unique, because we want duplicates
    substr(., 1, nchar(.) - 1)
  
  # add terminal parents by...
  # getting 
  terminalParents <- unknowns[union(which(duplicated(unknowns)), which(duplicated(unknowns, fromLast = T)))]
  names(terminalParents) <- terminalParents %>% 
    map(~ strsplit(.x, ".", fixed = T)[[1]] %>% last)
  terminalParents %<>% 
    map(~ 2 ^ -(1 + lengths(regmatches(.x, gregexpr(".", .x, fixed = T))))) %>%
    unlist(.) %>%
    tapply(., names(.), sum)
  
  # add unknown parents
  unknownParents <- setdiff(unknowns, unknowns[duplicated(unknowns)]) %>%
    map(~ 2 ^ -(1 + lengths(regmatches(.x, gregexpr(".", .x, fixed = T))))) %>%
    unlist(.) %>%
    sum()
  names(unknownParents) <- "Unknown"
  if(unknownParents == 0){unknownParents <- NULL}
  
  return(c(terminalParents, unknownParents) %>% signif(digits = 4) %>% sort(decreasing = T))
  
}