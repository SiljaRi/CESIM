# Functions (also) used by Phylogeny.R####

#' Util function that checks that CCFs stay 0 once a sample is 0 = NO Re-Appearance
#' 
#' @param clone CCFs of a clone
#' @return corrected CCFs of the clone (no reapperance)
correct_clone_disappear <- function(clone){
  gone <- FALSE
  pos <- FALSE
  for(i in 1:length(clone)){
    #Clone appeared (= prevalence > 0)
    if(clone[i]>0){
      pos <-TRUE
    } 
    #First Change to zero (Clone disappeared)
    if(pos == TRUE && clone[i] == 0){
      gone <- TRUE
    }
    #all following positions stay 0
    if(gone){
      clone[i] <- 0
    }
  }
  return(clone)
}


#' Util function that checks if clone CCFs are greater than Detection Minimum
#' 
#' @param clone CCFs of a clone
#' @param mindetec Minimum Detection value
#' @return correctet CCFs of a clone (only > 0 if greater than mindetec)
correct_clone_detection <-function(clone, mindetec){
  for(i in 1:length(clone)){
    if(clone[i]<mindetec){# CCF too small, in real data not detectable
      clone[i] = 0
    }
  }
  return(clone)
}

#' Util function that checks distance of clone to each other clone 
#' 
#' @param clone CCFs of a clone
#' @param table with CCFs of other clones (not containing CCF of param clone)
#' @param dist Minimun Distance value
#' @return TRUE if distance from clone to all other clones is > dist, else FALSE
check_clone_distance <- function(clone, table, dist){
  if(sum(clone)==0){return (FALSE)}
  if(is.vector(table)){
    if(length(clone)>1){#Several samples, only 1 former Clone so far
      if(!max(abs(table - clone))>=dist){
        return(list(FALSE, 1))
      }
    }else{#One Sample, table contains all clones
      for(i in 1:length(table)){
        if(!max(abs(table[i] - clone))>=dist){ # Difference smaller than dist
          return(list(FALSE, i))
        }
      }
    }
} else {# Table is still table (= At least 2 clones)
  for(i in 1:nrow(table)){
    if(!max(abs(table[i,] - clone))>=dist){ # Difference smaller than dist
      return(list(FALSE, i))
    }
  }
}
return(TRUE)
}

#' Util function that checks that sum(CCF) does not exceed 100 in total and in all clones
#' 
#' @param table with parents and CCFs of all clones (no Ancestor)
#' @return  list of T/F, and String: TRUE if sum of CCFs <= 100 and success mesasge, else FALSE with an Error message
check_frequencies <- function(table){
  # CCF of Clone i in sample t = Sum of all children(subclones) of clone i
  for(i in 1:nrow(table)){
    children <- which(table[,1]==i)
    if(length(children) == 0){next}
    for(j in 2:ncol(table)){
      sum <- 0
      for(k in 1:length(children)){
        sum = sum + table[children[k],j]
      }
      if((sum-table[i, j])>tol){ # Sum of chiildrens CCF are bigger than parent CCF - not possible in real data
        return(list(FALSE,paste(c("ERROR: Sum of all subclones can not be bigger than clone itself. 
                                  Problem in clone",i," with sublone(s)", children, "at sample",(j-1)),sep="",collapse =",")))
      }
    }
  }
  # Total CCF at sample t = sum of all independent clones at t
  ur <- which(table[,1]==0) # independent: parent = 0
  for(x in 2:ncol(table)){
    sum <- 0
    for(y in 1:length(ur)){
      sum = sum + table[ur[y],x]
    }
    if((sum - 100) > tol) {
      return(list(FALSE,paste("ERROR: There is more than 100% at timepoint", x-1, "with", sum, "\n")))
    }
  }
  return(list(TRUE, paste("Update successful")))
}

# Main Function for Proofing a valid Phylotable after update ####

#' Main function that makes several checks of phylogeny table after an update in the simulator
#' 
#' @param table phylogeny/clone table with parents, load and CCFs of all clones 
#' @param mindetec Minimum Detection value
#' @param mutationcount count of mutations
#' @param mindist Minimun Distance value
#' @param cnvs cnv table
#' @return list of T/F, and String: TRUE if all checks were positiv and possible warning message, else FALSE with an Error message
check_phylogeny <- function(table, mindetec, mutationcount, mindist, cnvs){
  Parent <- table[,1]
  Load <- table[,2]
  table <- table[,-(1:2)] # only CCFs
  warning <- NULL
  loadcheck <- check_load(Load, mutationcount, cnvs)
  if(!loadcheck[[1]]){
    return(list(FALSE, loadcheck[[2]]))
  }else{
    warning <- loadcheck[[2]]
  }
  if(is.vector(table)){# Only 1 Sample length = Anzahl Clone
    for(i in 1:length(table)){
      table[i] <- correct_clone_detection(table[i],mindetec)
      if(!check_parents(i, Parent[i], length(table))){
        return(list(FALSE, paste("Clone",i,"has invalid parent.")))
      }
      if(i == 1){ 
        next}
      simclone <-check_clone_distance(table[i],table[1:(i-1)],mindist)
      if(!simclone[[1]]){
        if(sum(table[i])==0){
          return(list(FALSE, paste("Clone",i,"lies under detection minimum or is zero in all samples.")))
        }
        return(list(FALSE, paste("Clone",i, "is too similar to clone", simclone[[2]],". (Minimum Difference is a frequency of 2 in at least one sample)")))
      }
    }
  }else{# Min 2 sample
    for(i in 1:nrow(table)){
      table[i,] <- correct_clone_detection(table[i,],mindetec)
      table[i,] <- correct_clone_disappear(table[i,])
      if(!check_parents(i, Parent[i], nrow(table))){
        return(list(FALSE, paste("Clone",i,"has invalid parent.")))
      }
      if(i == 1){ 
        next}
      simclone <-check_clone_distance(table[i,],table[1:(i-1),],mindist)
      if(!simclone[[1]]){
        if(sum(table[i,])==0){
          return(list(FALSE, paste("Clone",i,"lies under detection minimum or is zero in all samples.")))
        }
        return(list(FALSE, paste("Clone",i, "is too similar to clone", simclone[[2]],". (Minimum Difference is a frequency of", mindist,"in at least one sample)")))
      }
    }
  }
  if(!is.vector(table)){
    samplecheck <- check_sample(table)
    if(!samplecheck[[1]]){
      warning <- paste(warning, samplecheck[[2]])
    }
    if(!is.null(warning)){
      warning <- paste("WARNING: ", warning)
    }
  }
  freqs <- check_frequencies(cbind(Parent, table))
  if(freqs[[1]]){
    if(is.vector(table)){
      return(list(TRUE, paste(warning, freqs[[2]]), cbind(Parent, Mutations = Load, t1 = table)))
    }else{
      return(list(TRUE, paste(warning, freqs[[2]]), cbind(Parent, Mutations = Load, table)))
    }
  } else {
    return(list(FALSE, paste(freqs[[2]])))
  }
}

# Util Functions ####
#Hilfsfunktion, die Load überprüft
#' Util function that checks for a valid loadvector (each clone min. 1 mutation) and compares it to Number of Mutations in Basic Settings
#' 
#' @param loadvector vector with assignment of mutationcount to clones
#' @param mutationcount count of mutations
#' @param cnvs cnv table
#' @return list of T/F, and String: TRUE if loadvector is valid and possible warning, else FALSE and Error message
check_load <- function(loadvector, mutationcount, cnvs){
  sum <- 0
  for(i in 1:length(loadvector)){
    if(loadvector[i] <= 0){
      return(list(FALSE, paste("ERROR: Count of mutations has to be greater than 0 in clone", i, ".")))
    }
    sum <- sum + loadvector[i]
  }
  if(sum > mutationcount){
    return(list(TRUE, paste("Sum of mutations is greater than chosen mutation count.")))
  }else if(sum < mutationcount){
    return(list(TRUE, paste("Sum of mutations is smaller than chosen mutation count.")))
  }else {
    return(list(TRUE, NULL))
  }
}

#' Uil function that checks if a sample is not 0 for all clones
#' 
#' @param table phylogeny/clone table CCFs of all clones 
#' @return list of T/F and String: TRUE if there is at least one detection in each sample and no message, else FALSE and Error Message
check_sample <- function(table){
  for(i in 1:ncol(table)){
    if(sum(table[,i])==0){
      return(list(FALSE, paste("No detection in sample", i, ".")))
    }
  }
  return(list(TRUE, NULL))
}

#' Util function that checks for valid parent value
#' 
#' @param clonenumber index of clone to be tested
#' @param parent index of parent clone of clonenumber
#' @param count count of clones
#' @return TRUE if parent value is valid, else FALSE
check_parents <-function(clonenumber, parent, count){
  if(clonenumber == parent || parent > count){# parent invalid
    return(FALSE)
  }
  return(TRUE)
}

#' Util function that deletes a clone and corrects all other entries (ids, parents)
#' 
#' @param table phylogeny/clone table with parents, load and CCFs of all clones 
#' @param clonenumber index of clone to be deleted
#' @return corrected table without deleted clone
delete_clone <- function(table, clonenumber){
  table <- as.data.frame(table)
  children <- which(table[,1]==clonenumber)
  if(length(children) == 0){#Deleted Clone has no Subclone
    table <- table[-clonenumber,]
    for(i in 1:nrow(table)){
      if(table[i,1] >= clonenumber){# parent greater than index of deletec clone -> decrease
        table[i,1] <- (table[i,1]-1)
      }
    }
  } else {# Deleted Clone has sublcones
    for(i in 1:length(children)){
      currparent <- table[children[i],1]
      newparent <- table[currparent,1]
      table[children[i],1] <- newparent
    }
    table <- table[-clonenumber,]
    for(i in 1:nrow(table)){
      if(table[i,1] >= clonenumber){# parent greater than index of deletec clone -> decrease
        table[i,1] <- (table[i,1]-1)
      }
    }
  }
  rownames(table) <- 1:nrow(table)
  return(table)
}