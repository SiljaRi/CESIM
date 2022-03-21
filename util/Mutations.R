library(gtools)

#' Util function to get all possible clusters for an overlapping CNV
#'
#' @param freqvector vector with CCFs of CNV that should overlap
#' @param clones table with all clone CCFs
#' @param func function of CNV 
#' @param index clone containing the CNV
#' @return a vector with all possible clusters for the given overlapping CNV
get_clusters <- function(freqvector, clones, func, index){
  freqtable <- clones[,-(1:2)]/100 # clonefractions in 0-1
  clusterlist <- NULL
  
  #collect all sublones (of CNV) = children downwards
  parents <- clones[,1]
  children <- index
  for(i in 1:nrow(clones)){
    for(j in 1:length(children)){
      if(parents[i] == children[j]){
        children <- unique(c(children, i))
      }
    }
  }
  
  #collect all parents and grandparents (of CNV) upwards
  grandparents <- index
  latest <- index
  for(i in 1:nrow(clones)){
    if(latest == 0){
      break
    }
    latest <- parents[latest]
    grandparents <- c(grandparents, latest)
  }
  set <- 0:6
  setdiff(set, 0)
  grandparents <- setdiff(grandparents, 0)
  if(func == "CNV first"){#vector >= table
    for(i in 1:length(children)){
      for(j in 1:ncol(freqtable)){
        if(freqvector[j] < freqtable[children[i],j]){
          break
        } 
        if (j == ncol(freqtable)){
          clusterlist <- c(clusterlist, children[i])
        }
      }
    }
  } else if(func == "SNV first (active)" || func == "SNV first (passive)"){ #vector <= table
    for(i in 1:length(grandparents)){
      for(j in 1:ncol(freqtable)){
        if(freqvector[j] > freqtable[grandparents[i],j]){
          break
        } 
        if (j == ncol(freqtable)){
          clusterlist <- c(clusterlist, grandparents[i])
        }
      }
    }
  } else if(func == "Parallel"){ #neither grandparents nor children 
    all <- 1:nrow(clones)
    clusterlist <- setdiff(all,(union(children, grandparents)))
  }else{
    cat("Unknown Function Type:", func, "\n")
  }
  return(clusterlist)
}

#' Main function that creates all SNV data
#'
#' @param clones phylogeny table of simulator
#' @param coverage mean coverage value
#' @param purity purity vector
#' @param mindistance Minimum Distance value of simulator
#' @param cnvs CNV table of simulator
#' @return a list with the mutation datatable and a possible warning string
create_snvs <- function(clones, coverage, purity, mindistance, cnvs){
  clones <- as.data.frame(clones)
  parents <- as.vector(clones[,1])
  load <- as.vector(clones[,2])
  mutcount <- as.numeric(sum(load))
  cloncount <- as.numeric(nrow(clones))
  samcount <- as.numeric(ncol(clones)-2)
  coverage <- as.numeric(coverage)
  purity <- as.vector(purity)/100
  mindistance <- as.numeric(mindistance)
  
  ## Correct the loadvector by subtracting the count of cnvs that are already created as mutations
  if(!is.null(cnvs)){
    for(i in 1:nrow(cnvs)){
      clo <- cnvs[i,"in_Clone"]
      load[clo] <- load[clo] -1
    }
    cnvcount <- nrow(cnvs)
  } else {
    cnvcount <- 0
  }
  snvcount <- mutcount - cnvcount

  #create all columnnames
  refnames <- paste(rep("ref",times = samcount), as.character(1:samcount), sep = "_")
  varnames <- paste(rep("var",times = samcount), as.character(1:samcount), sep = "_")
  vafnames <- paste(rep("vaf",times = samcount), as.character(1:samcount), sep = "_")  
  ccfnames <- paste(rep("ccf",times = samcount), as.character(1:samcount), sep = "_")
  names <- vector()
  for(i in 1:samcount){
    names <- c(names, refnames[i], varnames[i], vafnames[i], ccfnames[i])
  }
  geno <- vector()
  inclone <- vector()
  incnv <- vector()
  
  #Start creating SNV
  if(snvcount > 0){
    mutcount <- 0
    reftable <- matrix(data = 0, nrow = 1, ncol = 4*samcount)
    id <- paste(rep("SNV", times = snvcount), as.character(1: snvcount), sep='_')
    chr <- vector()
    start <- vector()
    warning <- vector()
    end <- vector()
    # Build Overlap with CNV
    if(!is.null(cnvs)){
      for(i in 1:nrow(cnvs)){
        if(cnvs[i,"Overlap"]){
          count <- cnvs[i,"SNVs"]
          if(count < 1){
            next
          } else {# Overlap exists
            pchr <- cnvs[i, "Chr"]
            pstart <- cnvs[i, "Start"]
            pend <- cnvs[i, "End"]
            pFunction <- cnvs[i, "Function"]
            clist <- get_clusters(cnvs[i,-(1:9)], clones, pFunction, cnvs[i, "in_Clone"]) # get all possible clusters for overlapping CNV
            if(length(clist) == 0){# Overlap not possible
              warning <- c(warning, paste("WARNING: There is no solution to overlap CNV",i,"under the given Function",pFunction,"(no fitting ClusterCCFs). Overlap ignored."))
              next
            } 
            nosolution <- FALSE
            for(j in 1:count){# create overlapping SNV
              if(nosolution){
                break
              }
              if(length(clist)>1){
                clist <- sample(clist) #randomize order of clist
              }
              for(k in 1:length(clist)){
                if(load[clist[k]] > 0){
                  chr <- c(chr, pchr)
                  start <- c(start, sample(pstart:pend, 1))
                  clonenumber <- clist[k] # Assingment of SNV to a cluster from clist
                  inclone <- c(inclone, clonenumber)
                  incnv <-c(incnv, i)
                  ccfs <- as.numeric(clones[clonenumber,-(1:2)]/100) # Take according clusterfreqs
                  w <- x <- y <- z <- vector()
                  cnvvalue <- NULL
                  vafs <- vector()
                  # Set genotype and VAF values
                  switch(pFunction,
                         "CNV first" = {
                           if(cnvs[i, "Type"] == "Deletion"){
                             geno <- c(geno, "B")
                           } else if (cnvs[i, "Type"] == "Duplication"){
                             geno <- c(geno, "AAB")
                           } else if (cnvs[i, "Type"] == "LOH"){
                             geno <- c(geno, "AB(loh)")
                           }
                           cnvvalue = 1
                           for(m in 1:samcount){
                             ccfcnv <- cnvs[i, 9+m]
                             w[m] <- ccfs[m] #simple, no y present
                             x[m] <- ccfcnv-w[m] #is negative, if CCFSNV > CCFCNV
                             if(x[m] < 0 ){
                               cat("X NEGATIVE! (Should not happen)\n")
                             }
                             y[m] <- 0
                             z[m] <- 1 - ccfcnv
                           }
                         },
                         "SNV first (active)" = {
                           if(cnvs[i, "Type"] == "Deletion"){
                             geno <- c(geno, "AB(d)")
                             cnvvalue <- 0
                           } else if (cnvs[i, "Type"] == "Duplication"){
                             geno <- c(geno, "ABB/AB")
                             cnvvalue <- 2
                           } else if (cnvs[i, "Type"] == "LOH"){
                             geno <- c(geno, "BB/AB")
                             cnvvalue <- 2
                           }
                           for(m in 1:samcount){
                             ccfcnv <- cnvs[i, 9+m]
                             y[m] <- ccfs[m] - ccfcnv # Problem if CNV > SNV
                             x[m] <- 0
                             w[m] <- ccfcnv # simple no x present
                             z[m] <- 1 - ccfs[m]
                           }
                         },
                         "SNV first (passive)" = {
                           if(cnvs[i, "Type"] == "Deletion"){
                             geno <- c(geno, "B/AB")
                           } else if (cnvs[i, "Type"] == "Duplication"){
                             geno <- c(geno, "AAB/AB")
                           } else if (cnvs[i, "Type"] == "LOH"){
                             geno <- c(geno, "AB(loh)")
                             cnvvalue <- 0
                           }
                           cnvvalue = 1
                           for(m in 1:samcount){
                             ccfcnv <- cnvs[i, 9+m]
                             y[m] <- ccfs[m] - ccfcnv # Problem if CNV > SNV
                             x[m] <- 0
                             w[m] <- ccfcnv # simple no x present
                             z[m] <- 1 - ccfs[m]
                           }
                         },
                         "Parallel" = {
                           if (cnvs[i, "Type"] == "LOH"){
                             geno <- c(geno, "AB(loh)")
                           }else{
                             geno <- c(geno, "AB(p)")
                           }
                           cnvvalue = 1
                           for(m in 1:samcount){
                             ccfcnv <- cnvs[i, 9+m]
                             y[m] <- ccfs[m]
                             x[m] <- ccfcnv
                             w[m] <- 0
                             z[m] <- 1 - (ccfs[m]+ccfcnv)
                           }
                         }
                  )
                  switch(cnvs[i,"Type"],
                         "Deletion"= {
                           # Formula  Deletion
                           for(m in 1:samcount){
                             vafs[m] <- ((cnvvalue * w[m] + y[m])/(w[m] + x[m] + 2*y[m] + 2*z[m]))
                           }
                         },
                         "Duplication" = {
                           #Formula Duplication
                           for(m in 1:samcount){
                             vafs[m] <- (((cnvvalue * w[m]) + y[m])/(3*w[m] + 3*x[m] + 2*y[m] + 2*z[m]))
                           }
                         },
                         "LOH" = {
                           #Formula LOH
                           for(m in 1:samcount){
                             vafs[m] <- (((cnvvalue * w[m]) + y[m])/2)
                           }
                         }
                  )
                  mutcount <- mutcount + 1
                  load[clist[k]] <-  load[clist[k]]-1
                  reftable <- create_reads(samcount, coverage, vafs, mindistance, reftable, ccfs, purity)# Creating Reads
                  break
                } else if(k == length(clist)){
                  warning <- c(warning, paste("WARNING: There is no solution to overlap CNV",i,"under the given Function",pFunction,"with", count," mutations. (Mutationload in possible clones too few). Overlap (partially) ignored."))
                  nosolution <- TRUE
                  break
                }
              }
            }
          }
        }
      }
      end <- start
    }
    restsnvs <- snvcount - mutcount #remaining SNVs without Overlap
    if(restsnvs > 0){
      for(i in 1:restsnvs){
        #Create non-overlapping SNV
        repeat{
          no_overlap <- TRUE
          pchr <- sample.int(22, 1) # No sex chromosones!
          if(!is.null(cnvs)){
            for(j in 1:nrow(cnvs)){
              if(cnvs[j,"Chr"] == pchr){
                no_overlap <- FALSE
              }
            }
          }
          if(no_overlap){break}
        }
        chr <- c(chr, pchr)
        pstart <- sample.int(1000000, 1, prob = dnorm(1:1000000, mean = 5000, sd = 10000))
        start <- c(start, pstart)
        end <- c(end, pstart)
        geno <- c(geno, "AB")
      }
    }
    for(i in 1:cloncount){# iterate through all clones
      if(load[i] == 0){
        # no more mutations needed in this clone
        next
      }
      ccf <-as.numeric(clones[i,-(1:2)]/100)
      for(j in 1:load[i]){
        inclone <- c(inclone, i)
        incnv <- c(incnv, -1)
        vafs <- ccf/2
        reftable <- create_reads(samcount, coverage, vafs, mindistance, reftable, ccf, purity) # Create Reads for normal SNV
      }
    }
    # All SNVs are created. Connecting the data.
    mutationtable <- data.frame(ID = id, Chr = chr, Start = start, End = end, Geno = geno, Overlap_CNV = incnv, Clone = inclone )
    reftable <- as.data.frame(reftable)
    reftable <- reftable[-1,] # Cut empty util line
    colnames(reftable) <- names
    mutationtable <- cbind(mutationtable, reftable)
  } else {# No SNV only CNV
    mutationtable <- data.frame(ID = 0, Chr = 0, Start = 0, End = 0, Geno = NA, Overlap_CNV = NA, Clone = NA )
    reftable <- matrix(rep(0, times = samcount*8), nrow = 2)
    colnames(reftable) <- names
    mutationtable <- cbind(mutationtable, reftable)
    mutationtable <- mutationtable[-(1:2),] # provide empty table for showsnvs
  }
  return(list(mutationtable, warning))
}

#' Main function that creates read data for one SNV
#' 
#' @param samcount count of samples
#' @param coverage mean coverage value
#' @param vafs vector with calculated vafs (corrected for overlapping CNV) 
#' @param mindistance Minimum Distance value of simulator
#' @param reftable table with already created read data of previous SNVs
#' @param ccf CCFs of the SNV
#' @param purity purity vector
#' @return a list with the mutation datatable and a possible warning string
create_reads <- function(samcount, coverage, vafs, mindistance, reftable, ccf, purity){
  repeat{
    divisor <- sample(1:1000, samcount, replace = TRUE, prob = dlnorm(1:1000, meanlog = log(mean(coverage)), sdlog = 0.7, log = F)) #depth, with examined distribution
    vars <- round(vafs*divisor,0)
    refs <- divisor - vars
    # Conidering purity
    vars <- ceiling(vars * purity)
    check_difference<- FALSE
    # Check if reads fit to vaf inside minimum distance bound, so they are assigned to the right cluster
    if(mindistance > 0){
      for(k in 1:samcount){
        if(abs((vars[k]/(refs[k]+vars[k])) - vafs[k]) >= mindistance/2){
          print(abs((vars[k]/(refs[k]+vars[k])) - vafs[k]))
          check_difference <- TRUE
          cat("Differenz in Vars zu gering\n") # intern print
        }
      }
    }
    if(!check_difference){
        break}
  }
  row <- vector()
  for(k in 1:samcount){
    row <- c(row, as.numeric(refs[k]), as.numeric(vars[k]), as.numeric(round((vars[k]/(refs[k]+vars[k])),3)),
             as.numeric(round(ccf[k],3)))
  }
  reftable <- rbind(reftable, row) # Add new row to existing reftable
  return(reftable)
}


#' Main function that creates all CNV data
#'
#' @param cnvcheckbox the input of the explicit cnv checkbox in simulator
#' @param cnvcount count of CNVs
#' @param clones the CCFs of all clones
#' @return a table with CNV data
create_cnvs <- function(cnvcheckbox, cnvcount, clones){
  dup    <- "Duplication"  %in% cnvcheckbox
  del    <- "Deletion"     %in% cnvcheckbox
  loh    <- "LOH"     %in% cnvcheckbox
  count <- cnvcount
  ID <- paste(rep("CNV", times = count), as.character(1:count), sep='_')
  End <-vector()
  repeat{
    Chr <- sample.int(22, count, replace = TRUE) # No Sex Chromosomes accounted.
    Start <- sample.int(1000000, count, replace = TRUE, prob = dnorm(1:1000000, mean = 5000, sd = 10000))
    for(i in 1:count){
      End[i] <-  Start[i] + sample.int(10000, 1, replace = TRUE, prob = dnorm(1:10000, mean = 2000, sd = 10000))
    }
    if(check_overlap(data.frame(Chr, Start, End))[[1]]){break}
  }
  Type <- NULL
  if(dup){
    Type <- "Duplication"
  }
  if(del){
    Type <- c(Type, "Deletion")
  }
  if(loh){
    Type <- c(Type, "LOH")
  }
  # Collect |cnvcount| Types
  if(length(Type) < cnvcount){
    if(is.null(cnvcheckbox)){
      rest <- sample(c("Duplication", "Deletion", "LOH"), cnvcount-length(Type), replace = TRUE, prob = c(0.3,0.3,0.3))
      Type <- c(Type, rest)
    } else {
      rest <- sample(cnvcheckbox, cnvcount-length(Type), replace = TRUE)
      Type <- c(Type, rest)
    }
  } else {
    Type <- Type[1:cnvcount]
  }
  # Start creating CNV
  # CCF is taken in order from Clone 1 to Clone N
  # If more CNV than clones this scheme repeats 
  # = equal distribution of CNVs Clones
  in_Clone <- vector()
  if(nrow(clones)>=count){
    freqs <- clones[1:count,-(1:2)]
    in_Clone <- 1:count
  }else{
    freqs <- clones[1:nrow(clones), -(1:2)]
    in_Clone <- 1:nrow(clones)
    rest <- count - nrow(freqs)
    load <- clones[,2]-rep(1, times = nrow(clones))
    for(i in 1:rest){
      for(j in 1:length(load)){
        if(load[j]>0){
          freqs <- rbind(freqs, clones[j,-(1:2)])
          in_Clone <- c(in_Clone, j)
          load[j]<- load[j]-1
          break
        }
      }
    }
  }
  freqs <- freqs/100
  cnvs <- cbind(data.frame(ID, Chr, Start, End, Type, Overlap = FALSE, SNVs = rep(0, times = count)), Function = NA, in_Clone, freqs)
  return(cnvs)
}


#' Main function that transforms CNV data into format for mutation table
#'
#' @param cnvs cnv table (including type, function etc which is not shown in mutation table)
#' @return cnv data in format of mutation table
create_showcnvs <- function(cnvs){
  showcnvs <- NULL
  for(i in 1:nrow(cnvs)){
    row <- vector()
    for(j in 1:(ncol(cnvs)-9)){#jedes Sample
      row <- cbind(row, NA, NA, NA, (cnvs[i,(9+j)]))## Reftable CNVS ohne Reads und VAF (3*NA) nur freqs
    }
    completerow <- cbind(cnvs[i,1:4], NA, NA, cnvs[i,9],row)
    showcnvs <- rbind(showcnvs, completerow)
  }
  return(showcnvs)
}
