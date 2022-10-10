# Calculating completeness metrics of Triton Database

rm(list=ls())

library(divDyn)

group_names <- c("forams")

bin_sizes <- c(100000, 200000, 500000, 1000000)

bin_sizes_strings<-c("100000", "200000", "500000", "1000000")


# prepare empty dataframe to be filled
completeness_data <- data.frame(matrix(nrow = length(bin_sizes), ncol = 6))

colnames(completeness_data) <- c("group", "binsize", "bin_avg_proportion_part_timer", "all_bin_avg_three_timer","three_timer_SD", "whole_dataset_SCM")


fi<-1

dev_df_count <- 0


  ori_dat <- triton.pres #This object comes from Triton downloadable file
  
  
  #########################################################################
  # 2. taxonomic filtering
  
  
  # Let's make sure every occurrence has a species 
  sum(is.na(ori_dat$species))
  # Looking good (no NA's)
  
  # Omit non-informative species entries
  ori_dat <- ori_dat[ori_dat$species!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$species)==F, ]
  
  # Omit non-informative paleo-coordinate entries
  ori_dat <- ori_dat[ori_dat$pal.lat!="", ]
  ori_dat <- ori_dat[ori_dat$pal.long!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$pal.lat)==F, ]
  ori_dat <- ori_dat[is.na(ori_dat$pal.long)==F, ]
  
  
  
  # Taxonomic subset 
  
  # For now, keep everything
  
  
  #########################################
  
  # I will make a genus_species column, in case a species name is shared between two genera:
  
  ori_dat$species <- sub(" ", "_", ori_dat$species)
  
  #head(ori_dat$gen_sp)
  
  #set this column to whatever we want to do analysis on (genus or species), so we don't have to change the rest of the code
  ori_dat$taxon <- ori_dat$species
  #ori_dat$taxon <- ori_dat$Genus
  
  head(ori_dat$taxon)
  
  # Exclude occurrences with unassigned ages 
  # I suppose this is redundant, since we drop the most recent bin from model analysis anyways
  ori_dat <- ori_dat[is.na(ori_dat$age) == F,]
  
  
  # Assign each occurrence to a bin
  
  
  for (size in 1:length(bin_sizes)){
    
    dev_df_count <- dev_df_count + 1
    
    print(dev_df_count)
    
    # select bin size
    binsize <- bin_sizes[size]
    
    nb=ceiling(max(ori_dat$age)/(binsize/1000000))
    # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
    
    
    labs <- c(1:nb)
    
    labs <- as.character(labs)
    
    ori_dat$bin <- as.numeric(cut(ori_dat$age, breaks=(nb), labels = labs))
    
    # divDyn automatically calculates several metrics of interest
    out <- divDyn(ori_dat, tax = "taxon", bin = "bin")
    
    completeness_data[dev_df_count, 1] <- group_names[fi]
    
    completeness_data[dev_df_count, 2] <- bin_sizes[size]
    
    # part-timers divided by all range-through taxa, not sure if this is accepted in the literature, but I'm experimenting with this one
    completeness_data[dev_df_count, 3] <- 1-mean(out$tPart/out$tThrough, na.rm = T) 
    
    # three timer completeness (calcuated for each bin), averaged accross all bins for a organism type/bin size combination
    completeness_data[dev_df_count, 4] <- mean(out$samp3t, na.rm=T)
    
    # standard deviation of those averages for all bins
    completeness_data[dev_df_count, 5] <- sd(out$samp3t, na.rm=T)
    
    
    
    # simple completeness metric (SCM)
    # Ratio of observed occurrences (in whole data set) to total inferred occurrences (in whole data set),
    # where an occurence is a given taxon being present in a given time bin (i.e. multiple fossil occurences of the 
    # same taxon in the same time bin only counts as one occurrence)
    
    # see Fara and Benton (2000)
    
    completeness <- data.frame(matrix(ncol = 2, nrow = (length(unique(ori_dat$taxon))+9)))
    
    colnames(completeness) <- c("all_inferred_occurences", "observed_occurences")
    
    for (i in 1:length(unique(ori_dat$taxon))){
      
      tmp <- (unique(ori_dat$taxon))[i] # for each taxon
      
      tmpbin <- sort(unique(ori_dat$bin[ori_dat$taxon==tmp])) #get list of all time bins in which that taxon occurs
      
      binmin <- min(unique(ori_dat$bin[ori_dat$taxon==tmp]))
      
      binmax <- max(unique(ori_dat$bin[ori_dat$taxon==tmp]))
      
      rng <- binmax - binmin + 1
      
      completeness[i,1] <- rng # total inferred fossil occurrences for each taxon
      
      completeness[i,2] <- length(tmpbin) # all recorded time bins for each taxon (observed occurrences) 
      
    }
    
    completeness_data[dev_df_count, 6] <- (sum(completeness$observed_occurences, na.rm=T))/(sum(completeness$all_inferred_occurences, na.rm=T))
    
    
    
    
    
  }
  
} 
