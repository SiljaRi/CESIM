#' Util function that checks if any two CNV overlap
#' 
#' @param cnvtable table with Chr, Start and End of all cnvs
#' @return list of T/F and string: TRUE if no overlap was detected and empty message, else FALSE and error message
check_overlap <- function(cnvtable){
  if(nrow(cnvtable) == 1){ # only one cnv, no overlap possbile
    return(list(TRUE, NULL))
  }
  for(i in 1:(nrow(cnvtable)-1)){
    for(j in (i+1):nrow(cnvtable)){
      if(cnvtable[i,1] == cnvtable[j,1]){# same Chr
        if(cnvtable[i,2] == cnvtable[j,2] || cnvtable[i,3] == cnvtable[j,3] || cnvtable[i,2] == cnvtable[j,3] || cnvtable[i,3] == cnvtable[j,2]){
          return(list(FALSE, paste("CNVs are not allowed to overlap. Problem with CNVs",i,",",j,".")))
        } else if(cnvtable[i,2] < cnvtable[j,2]){# i is smaller, has to end before j
          if(cnvtable[i,3] > cnvtable[j,2]){
            return(list(FALSE, paste("CNVs are not allowed to overlap. Problem with CNVs",i,",",j,".")))
          }
        }else if(cnvtable[i,2] > cnvtable[j,2]){ # i is bigger
          if(cnvtable[j,3] > cnvtable[i,2]){
            return(list(FALSE, paste("CNVs are not allowed to overlap. Problem with CNVs",i,",",j,".")))
          }
        } 
      }
    }
  }
  return(list(TRUE, NULL))
}

#' #' Util function that checks for enough SNV to be overlapped by CNVs
#' 
#' @param cnvtable table of all cnvs (cnv table of simulator)
#' @param snvcount count of SNVs
#' @return list of T/F and string: TRUE if enough SNV for overlap exist and possible warning, else FALSE and error message
check_snvcount <- function(cnvtable, snvcount){
  sum <- 0
  warning <- NULL
  for(i in 1:nrow(cnvtable)){
    if(cnvtable[i, "Overlap"]){
      if(cnvtable[i, "SNVs"] == 0){
        warning <- paste("WARNING: Overlap with 0 SNVs will be ignored.")
      }
      sum <- sum + cnvtable[i, "SNVs"]
    }
  }
  if(sum > snvcount){
    return(list(FALSE, paste("There are not enough SNVs to be overlapped.")))
  } else {
    return(list(TRUE, warning))
  }
}


#' #' Util function that checks for valid Start and End values of CNVs
#' 
#' @param cnvtable table of all cnvs (cnv table of simulator)
#' @return list of T/F and string: TRUE if positon is valid and empty message, else FALSE and error message
check_position <- function(cnvtable){
  for(i in 1:nrow(cnvtable)){
    if(cnvtable[i, "Start"] > cnvtable[i, "End"]){
      return(list(FALSE, paste("ERROR: End before Start at CNV", i,".")))
    }
  }
  return(list(TRUE, NULL))
}


#' #' Util function that checks for valid in_Clone values
#' 
#' @param cnvtable table of all cnvs
#' @param loadvector vector with assignments of mutationcounts to clones
#' @return list of T/F and string: TRUE if in_Clone value is valid and empty message, else FALSE and error message
check_inClone <- function(cnvtable, loadvector){
  for(i in 1:nrow(cnvtable)){
    if(cnvtable[i, "in_Clone"] == 0){
      return(list(FALSE, paste("ERROR: in_Clone can not be 0 in CNV", i, ".")))
    } else if(cnvtable[i, "in_Clone"] > length(loadvector)){
      return(list(FALSE, paste("ERROR: in_Clone value of CNV", i, " is not valid.")))
    } else {
      clone <- cnvtable[i, "in_Clone"]
      loadvector[clone] <- loadvector[clone] - 1
      if(loadvector[clone] < 0){
        return(list(FALSE, paste("ERROR: Not enough mutations in Clone", clone, "for CNV", i, ".")))
      }
    }
  }
  return(list(TRUE, NULL))
}

#' Util function that checks for valid function of any CNV
#' 
#' @param cnvtable table with all cnvs
#' @return list of T/F and string: TRUE if valid functions are given and empty message, else FALSE and error message
check_function <- function(cnvtable){
  for(i in 1:nrow(cnvtable)){
    if(cnvtable[i, "Overlap"]){
      if(is.na(cnvtable[i, "Function"])){
        return(list(FALSE, paste("ERROR: Please choose a function for overlapping CNV", i, ".")))
      }
    }
  }
  return(list(TRUE, NULL))
}


#' Main Check Function, checks CNV table for validity.
#' @return list of TRUE/NULL or FALSE/Errormessage, if an error occured while checking.
#' 
#' #' Main function that makes several checks on CNV table
#' 
#' @param cnvtable table of all cnvs
#' @param loadvector vector with assignments of mutationcounts to clones
#' @param snvcount count of SNVs
#' @param clonefreqs CCFs of all clones
#' @return list of T/F,string,dataframe: TRUE all checks were positive with possible warning and updated cnv table, else FALSE and error message (no table)
check_cnvs <- function(cnvtable, loadvector, snvcount, clonefreqs){
  if(nrow(cnvtable) == 0){
    return(list(TRUE, NULL, cnvtable))
  }
  pos <- check_position(cnvtable) # Check position
  if(!pos[[1]]){
    return(list(FALSE, pos[[2]]))
  } 
  lap <- check_overlap(cnvtable[,-1]) # Check overlap
  if(!lap[[1]]){
    return(list(FALSE, lap[[2]]))
  } 
  inclo <- check_inClone(cnvtable, loadvector) # Check in_Clone
  if(!inclo[[1]]){
    return(list(FALSE, inclo[[2]]))
  } 
  snvs <- check_snvcount(cnvtable, snvcount) # Check SNVcount
  if(!snvs[[1]]){
    return(list(FALSE, snvs[[2]]))
  } 
  func <- check_function(cnvtable) # Check Function
  if(!func[[1]]){
    return(list(FALSE, func[[2]]))
  }
  # Update CCFS for each CNV regarding in_Clone value
  for(i in 1:nrow(cnvtable)){
    if(is.vector(clonefreqs)){
       cnvtable[i,10] <- clonefreqs[cnvtable[i, "in_Clone"]]/100
    } else {
      for(j in 1:ncol(clonefreqs)){
        cnvtable[i,j+9] <- clonefreqs[cnvtable[i, "in_Clone"], j]/100
      }
    }
  }
  return(list(TRUE, snvs[[2]], cnvtable))
}