tol <- 1e-5 #threshold for

#' Main function that creates the phylogeny (table) 
#' 
#' @description This function is the main part of simulation. 
#'
#' @params all params (except counter) are inputs of the Basic Settings of the simulator
#' @param mutations count of mutations
#' @param clones count of clones
#' @param samples count of samples
#' @param evol kinf of evolution
#' @param mindetec Minimum Detection value
#' @param counter util parameter for counting the attempts of building a phylogeny
#' @param mindist Minimum distance value
#' @return a table with the created phylogeny (clone table)
create_phylogeny <- function(mutations, clones, samples, evol, mindetec, counter, mindist){
  mutations <- as.numeric(mutations)
  clones <- as.numeric(clones)
  samples <- as.numeric(samples)
  evol <-as.character(evol)
  counter <- as.numeric(counter)
  if(counter == 700){# 700 function calls on stack, abort simulation -> No phylogeny found
    return(NULL)
  }
  # Start building the first Clone (CCFs for all samples)
  repeat{
    Init_Clone <- round(runif(samples, min = mindetec, max = 100), 2)
    Init_Clone <- correct_clone_detection(Init_Clone, mindetec) # Check CCFs are > mindetec
    Init_Clone <- correct_clone_disappear(Init_Clone) # Check that Clone does not reappear
    if(sum(Init_Clone)!=0){#If Min Detection is high, correction can set clone to 0 in all samples.
      break
    }
  }
  #Table of Fractions (Ancestral (100), Clone (Init) , Restfraction (100-Init))
  Sampletable <- matrix(data = c(rep(x = 100, times = samples), Init_Clone, 100-Init_Clone),
                        nrow = 3, ncol = samples, byrow = TRUE)
  Sampletable <- check_minus(Sampletable)
  names <- c("Ancest", "Clone_1", "R_1")
  numbers <- c(0, 1, 0)
  parents <- c(0, 0, 0)
  children <- c(NA, NA, NA)
  splittable <- c(FALSE, check_split(Sampletable[2,]), check_split(Sampletable[3,])) # Vector giving all possible clones that can be split into new ones
  k <- 1
  # Create the N-1 remaining clones
  while(k < clones){
    if(sum(splittable)==0){
      # Not happened yet
      print("No more space to create new clones.\n")
      return(NULL)
    }
    count <- 0
    repeat{ 
      # One new clone
      splits <- which(splittable)
      # Decide on kind of evolution which clones are splitted
      if(evol == "linear"){
        splits <- subset(splits, splits %% 2 == 0)
        if(length(splits)>0){
          splitindex <- get_splitindex(splits)
        }else {
          cat("No linear separation possible. Start new simulation.\n")
          # Since values are created randomly, this run was bad random values. Just try a new run
          # Rerun function, increase counter
          return(create_phylogeny(mutations, clones, samples, evol, mindetec, (counter+1), mindist))
        }
      }else if(evol == "parallel-dependent"){# split cant be in 3
        par_splits <- setdiff(subset(splits, splits %% 2 == 1),3)
        splits <- c(par_splits, subset(splits, splits %% 2 == 0))
        if(length(splits) == 0 || (k == 2 && length(par_splits) == 0)){
          return(create_phylogeny(mutations, clones, samples, evol, mindetec, (counter+1),mindist))
        }
        if(k == 2){# 3. Clone is the first one possible
          splitindex <- get_splitindex(par_splits)
        }else{
          splitindex <- get_splitindex(splits)
        }
      }else if(evol == "independent"){# Independent = parent 0, else inear 
        # all uneven splittables, which have parent[x] == 0 
        indep_splits <- intersect(subset(splits, splits %% 2 == 1), which(parents == 0))
        splits <- c(indep_splits, subset(splits, splits %% 2 == 0))
        if(length (splits)==0 || (k == 1 && length(indep_splits) == 0)){
          return(create_phylogeny(mutations, clones, samples, evol, mindetec, (counter+1),mindist))
        }
        if(k == 1){# enforce independent Clone
          splitindex <- get_splitindex(indep_splits)
        }else{# every clone afterwards randomly linear or independent
          splitindex <- get_splitindex(splits)
        }
      }else { # random evolution choosen
        splitindex <- sample(x = which(splittable), size = 1)
      }
      New_Clone <- sapply(Sampletable[splitindex,], function(c){
        round(runif(1, min = 0, max = c), 2)}) # create clone CCFs
      New_Clone <- correct_clone_detection(New_Clone, mindetec)
      New_Clone <- correct_clone_disappear(New_Clone)
      # Test new clone
      Freqtable <- cbind(parents, Sampletable) #Sampletable with Parents for testing
      if(splitindex %% 2 == 0){ # EVEN: Clone descented from another clone
        row <- c(splitindex/2, New_Clone, splitindex/2, Sampletable[splitindex,]-New_Clone)
      }else{ # UNEVEN: Clone decended from restfraction
        row <- c(parents[splitindex], New_Clone, parents[splitindex], Sampletable[splitindex,]-New_Clone)
      }
      Freqtable <- rbind(Freqtable[-1,], matrix(row, nrow = 2, byrow = TRUE))
      Freqtable <- check_minus(Freqtable)
      Freqtable <- Freqtable[c(TRUE,FALSE),]
      if(check_clone_distance(New_Clone, Sampletable[c(FALSE, TRUE),], mindist)[[1]] && check_frequencies(Freqtable)[[1]]){break}
      count <- count + 1
      if(count > 100){
        # After 100 failed attempts to build new Clone, start new simulation
        return(create_phylogeny(mutations, clones, samples, evol, mindetec, (counter+1),mindist))
      }
    }
    # New Cone built successfully
    splittable[splitindex] <- FALSE
    Sampletable <- rbind(Sampletable, New_Clone, Sampletable[splitindex,]-New_Clone)
    Sampletable <- check_minus(Sampletable)
    splittable <- c(splittable, check_split(New_Clone), 
                    check_split(Sampletable[splitindex,]-New_Clone))
    numbers <- c(numbers, k+1, k+1)
    if(splitindex %% 2 == 0){ # EVEN: Clone descented from another clone
      parents <- c(parents, splitindex/2, splitindex/2)
      children[splitindex] <- k+1 
    }else{ # UNEVEN: Clone decended from restfraction
      parents <- c(parents, parents[splitindex], parents[splitindex])
    }
    children <-c(children, NA, NA)
    names <- c(names, paste0("Clone_",k+1), paste0("R_",k+1))
    k<-k+1
    
  }
  # All clones created, build output
  loadvector <- as.vector(get_loadvector(mutations, clones))
  resi <- data.frame(c(0, rep(1:clones, each = 2)), parents, loadvector, Sampletable)
  colnames(resi) <- c("Number","Parent", "Mutations", paste0("t", 1:samples))
  resi <- rbind(resi[1,],resi[ !c(TRUE,FALSE),])
  rownames(resi) <- 0:clones
  return (resi)
}


#' Util function that checks if clone can be split up (CCF > 5)
#' 
#' @param part CCF vector of clone
#' @return T if at least one CCF in one sample is greater 5, else F
check_split <- function(part){
  if(sum(part>=5.0)>=1){
    return(TRUE)
  }
  return(FALSE)
}


#' Util function that chooses random splittable clone to be split up
#' 
#' @param splitvector vector with all clones that can be split up
#' @return index of clone which will be split
get_splitindex <- function(splitvector){
  if(length(splitvector)==0){
    return(NULL)
    }
  if(length(splitvector)==1){
    return(splitvector)
  }else{
    return(sample(x = splitvector, size = 1))
  }
}


#' Util function that removes negative numbers which might be created by float-calculations
#' 
#' @param table table with possibly negative CCFs
#' @return table without negative CCFs (corrected to 0)
check_minus <- function(table){
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      if(table[i,j] < 0){
        table[i,j] = 0.00
      }
    }
  }
  return(table)
}

#' Util function that creates a randomized loadvector = assignment of mutationcount to different clones
#' 
#' @param mutations count of mutations
#' @param clones count of clones
#' @return randomized assingment of mutationcount to clones
get_loadvector <- function(mutations, clones){
  to_assign <- (mutations - clones)
  loadvector <- rep(1, times = clones)
  while(to_assign>0){
    index <- sample(1:clones, 1)
    load <- sample(1:to_assign, 1, prob = dnorm(1:to_assign, 3, 3))
    loadvector[index] <- (loadvector[index]+load)
    to_assign <- (to_assign - load)
  }
  loadvector <- c(0, rep(loadvector, each = 2))
  return(loadvector)
}

#Freqtable = Sampletable + neuer Clone für Checks
####| parents |   V2  |   V3  |  V4
#--------------------------------------
# 1 |    0    | 97.02 | 48.32 | 4.53      Clone
# 2 |    0    | 3.98  | 51.67 | 96.47     Rest
# 3 |    1    | 97.02 | 48.32 | 4.53      Clone 
# 4 |    1    |   0   |   0   |   0       Rest



#create_phylogeny(5,3,"independent",1, 0)

#



# TIMESCAPE INPUT ERZEUGEN
# phyl <- tabl[2:(clones+1),c(2,1)]
# colnames(phyl) <- c("source","target")
# 
# prev <- tabl[c(1:(clones+1)),c(3:(3+samples-1))]
# 
# prevt <- transpose(prev)
# 
# test <- data.frame(
#   timepoint = rep(1:samples, times = (clones+1)),
#   clone_id = rep(0:clones, each = samples),
#   clonal_prev = unlist(prevt)
# )
# 
# clone_colours <- data.frame(clone_id = c(0:clones), 
#                             colour = c("F6F6F6", grDevices::rainbow(clones)))
# 
# 
# timescape(clonal_prev = test, tree_edges = phyl, 
#           clone_colours = clone_colours,
#           height=350, alpha=15, 
#           genotype_position="space")
# 


