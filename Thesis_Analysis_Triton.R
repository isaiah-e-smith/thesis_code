# Isaiah E. Smith 
# Analysis of Triton Database data

rm(list=ls())

getwd()
setwd("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/")

library(gtools)
library(modEvA)
library(OneR)
library(tibble)
library(icosa)
library(divDyn)
library(chronosphere)
library(raster)
library(rgplates)



group_names <- c("forams")

bin_sizes <- c(100000, 200000, 500000, 1000000)

bin_sizes_strings<-c("100000", "200000", "500000", "1000000")

# create icosahedral grid, larger number = smaller grid cells
gr <- hexagrid(20, sp=TRUE)

# select how big of a bin gap change in occupancy should be calculated over (1 is default of consecutive bins)
intrvl <- 1

# create empty data frames to fill with percent deviance values 
deviance_explained_delta1 <- data.frame(matrix(ncol = 3, nrow = (length(bin_sizes))))
colnames(deviance_explained_delta1) <- c("group", "bin_size", "value")

a_deviance_explained_delta1 <- data.frame(matrix(ncol = 3, nrow = (length(bin_sizes))))
colnames(a_deviance_explained_delta1) <- c("group", "bin_size", "value")



#########################################################################
# 1. Load and explore the data

ori_dat <- triton.pres

#maximum age in data set
max(ori_dat$age, na.rm=T)

#number of records in data set
length(ori_dat$species)

#number of unique species in data set
length(unique(ori_dat$species))



#########################################################################
# 2. taxonomic filtering


# Let's make sure every occurrence has a Genus now
sum(is.na(ori_dat$species))
# Looking good (no NA's)

# Omit non-informative Genus entries
ori_dat <- ori_dat[ori_dat$species!="", ]
ori_dat <- ori_dat[is.na(ori_dat$species)==F, ]

# Omit non-informative paleo-coordinate entries
ori_dat <- ori_dat[ori_dat$pal.lat!="", ]
ori_dat <- ori_dat[ori_dat$pal.long!="", ]
ori_dat <- ori_dat[is.na(ori_dat$pal.lat)==F, ]
ori_dat <- ori_dat[is.na(ori_dat$pal.long)==F, ]



#########################################

# Make a genus_species column, in case a species name is shared between two genera:


head(ori_dat$species)

length(ori_dat$species)

length(unique(ori_dat$species))

# make new column called "taxon"
# set this column to whatever we want to do analysis on (genus or genus_species column)
ori_dat$taxon <- ori_dat$species

head(ori_dat$taxon)

# Exclude occurrences with unassigned ages 
ori_dat <- subset(ori_dat, is.na(age)==F)

length(ori_dat$species)

length(unique(ori_dat$species))


# Assign each occurrence to a bin

# counter to fill output dataframe
dev_df_count <- 0

# loop through our selected bin sizes
for (size in 1:length(bin_sizes)){
  
  dev_df_count <- dev_df_count + 1
  
  binsize <- bin_sizes[size]
  
  nb=ceiling(max(ori_dat$age)/(binsize/1000000))
  # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
  
  # create labels
  labs <- c(1:nb)
  
  labs <- as.character(labs)
  
  # use those labels to bin data
  ori_dat$bin <- as.numeric(cut(ori_dat$age, breaks=(nb), labels = labs))
  
  
  # As of now, the more recent bins have smaller numbers; 
  # i.e. as we progress forward through time, the bins get smaller (the present is bin 1)
  
  
  # Identify two-timers
  # Which taxa are recorded in at least two bins
  
  t <- unique(subset(ori_dat, select=c(taxon, bin)))
  t.2 <- tapply(t$bin, t$taxon, length)
  length(t.2)
  t.2 <- t.2[t.2 > 1]
  length(t.2)
  t2.names <- names(t.2)
  t2.mat <- data.frame(subset(ori_dat, taxon %in% t2.names)) # analysis matrix (dataframe)
  t2.inv <- unique(subset(t2.mat, select=c(taxon, bin)))
  
  
  # stratigraphic ranges of taxa
  orig <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, max)
  ext <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, min)
  
  
  # locate paleolatitude and paleolongitude coordinates on grid
  t2.mat$cells <- locate(gr, t2.mat[, c("pal.long", "pal.lat")])
  
  
  # make empty data frame
  un_factor <- data.frame(matrix(ncol = 6, nrow = nrow(t2.mat)))
  
  
  
  # For calculating raw in-bin proportion (without considering cells sampled in previous bin)
  # and with "bridging gap" for missing intervals
  
  count <- 0
  
  
  for (i in 1:length(unique(t2.mat$taxon))) { #for each specific taxon
    
    # specific taxon
    tx <- unique(t2.mat$taxon)[i]
    
    # creat list of bins in which that taxon occurs
    binlist <- sort(unique(t2.mat$bin[t2.mat$taxon==tx]))
    
    
    #for every occupied bin for each taxon
    for (j in 1:length(binlist)) {
      
      count <- count + 1
      
      if((count)%%100 == 0){
        print(count)
      }
      
      #specific bin from bin list
      bin <- binlist[j]
      
      prev_bin <- binlist[j+intrvl]
      
      # all geographic cells occupied (by any species) during this given bin
      a <- unique(t2.mat$cells[t2.mat$bin == bin])
      
      # all geographic cells occupied (by any species) during last occupied bin
      b <- unique(t2.mat$cells[t2.mat$bin == prev_bin])
      
      
      #now we can work with unique taxon/bin pairings
      
      # what unique cells are occupied by this taxon at this time bin
      un<-unique(t2.mat$cells[t2.mat$taxon== tx & t2.mat$bin == bin])
      
      # what unique cells are occupied by this taxon at the last occupied time bin (i.e. skip empty bins and "bridge gap")
      un_before <- unique(t2.mat$cells[t2.mat$taxon == tx & t2.mat$bin == (prev_bin)])
      
      # find raw proportion for this taxon in this bin
      rp <- length(un)/length(a)
      
      # find raw proportion for this taxon in previously occupied bin (NA if there is no "before" bin, i.e. if species just emerged)
      if (is.na(un_before)==T | is.na(b)==T | length(un_before)==0 | length(b)==0){
        rpb <- NA
      }else{
        rpb <- length(un_before)/length(b)
      }
      
      
      # proportion of present occupancy to previous occupancy
      prop <- (rp)/(rpb)
      
      # only keep if not NA or 0
      if (is.infinite(prop)==T) {
        ans <- NA
      } else if (is.na(prop) == T){
        ans <- NA
      } else if (prop==0){
        ans <- NA 
      } else {
        ans <- prop
      }
      
      
      # take natural log of the proportion value
      lans <- log(ans)
      
      # save to data frame
      un_factor[count, 1] <- bin #bin
      
      un_factor[count, 2] <- tx #taxon
      
      un_factor[count, 3] <- lans #log value
      
      un_factor[count, 4] <- rp #raw proportional occupancy
      
      un_factor[count, 6] <- length(un) #number of sampled cells for that taxon in that bin
      
    }
    
  }
  
  
  # name columns
  colnames(un_factor) <- c("bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled")
  
  # trim to only include rows that were filled
  un_factor <- un_factor[1:count,]
  
  
  
  # Fill Extinction column (1 = goes extinct or 0 = does not go extinct)
  
  # start by filling the Extinction column with all 0's (default: not going extinct in given bin)
  un_factor$extinction <- 0
  
  # Now, for last occurrences (extinctions) of each taxon, assign Extinction value of 1
  count <- 0
  
  for (i in 1:nrow(un_factor)) {
    
    count <- count+1
    
    if((count)%%10 == 0){
      print(count)
    }
    
    # makes list of occurrences for each taxon 
    # (this is currently redundant because as it iterates down the list, taxon will come up multiple times)
    # But it gets the job done
    spc <- un_factor[un_factor$taxon == un_factor$taxon[i],]
    
    # Redefines the Extinction variable as 1 (going extinct) for the first numerical bin occurrence (last temporal occurrence) of the taxon
    # i.e. where on the list the last occurrence of the taxon is
    un_factor$extinction[un_factor$bin == spc$bin[1] & un_factor$taxon == spc$taxon[1]] <- 1
    
    
  }
  
  
  
  #########################################
  # save the processed data to unique variable
  
  
  # Remove last bin (in this case, bin 1, 0.1 MYA-recent), because everything will automatically have an Extinction value of 1, 
  # since there are no further bins, even if they occur in the present and are not extinct
  un_factor_trim <- un_factor[!(un_factor$bin == min(un_factor$bin)),]
  
  # remove entries that have NA's for delta column (first occurences of species)
  # for longer intervals of occupancy change (i+2, i+3, etc... the number of NA entries for each species increases)
  un_factor_trim<-un_factor_trim[is.na(un_factor_trim$delta_1) == F,]
  
  # Save processed data to variable
  assign(paste("un_factor_trim", "NSB", group_names[fi], bin_sizes_strings[size], sep = "_"), (un_factor_trim))
  
  # make model
  model1 <- glm(extinction ~ raw_prop_sampled * delta_1, 
                family = binomial(link ="logit"), data = un_factor_trim)
  
  # step model
  step_model1 <- step(model1)
  
  # Save model output
  assign(paste("step_mod", "NSB", group_names[fi], bin_sizes_strings[size], sep = "_"), (step_model1))
  
  # Save ANOVA information
  af<-anova(step_model1)
  
  afss <- af$Deviance # deviance reduction attributed to each term
  
  af_prop<-(cbind(af,PrpExp=afss/sum(afss, na.rm=T)))
  
  # deviance reduction attributed to Delta_1
  PrpExp_D1 <- af_prop$PrpExp[rownames(af_prop) == "delta_1"]
  
  
  # Save this information to data frame
  deviance_explained_delta1$group[dev_df_count] <- paste(group_names[fi])
  deviance_explained_delta1$bin_size[dev_df_count] <- paste(bin_sizes[size])
  deviance_explained_delta1$value[dev_df_count] <- paste(PrpExp_D1)
  
  # same model but with different order of terms
  
  model1a <- glm(extinction ~ delta_1 * raw_prop_sampled, 
                 family = binomial(link ="logit"), data = un_factor_trim)
  
  step_model1a <- step(model1a)
  
  # Save ANOVA information
  afa<-anova(step_model1a)
  
  afssa <- afa$Deviance # deviance reduction attributed to each term
  
  af_propa<-(cbind(afa,PrpExp=afssa/sum(afssa, na.rm=T)))
  
  # deviance reduction attributed to delta_1
  PrpExp_D1a <- af_propa$PrpExp[rownames(af_propa) == "delta_1"]
  
  
  # Save this information to data frame
  a_deviance_explained_delta1$group[dev_df_count] <- paste(group_names[fi])
  a_deviance_explained_delta1$bin_size[dev_df_count] <- paste(bin_sizes[size])
  a_deviance_explained_delta1$value[dev_df_count] <- paste(PrpExp_D1a)
  
}


deviance_explained_delta1$value <- as.numeric(deviance_explained_delta1$value)
a_deviance_explained_delta1$value <- as.numeric(deviance_explained_delta1$value)


# add column to identify which formula order it came frome
deviance_explained_delta1$order <- "order_b"
a_deviance_explained_delta1$order[dev_df_count+1:]<- "order_a"


###############################################


# the 16 models (four NSB data sets and four bin sizes each)

all_calc_nano <- c("step_mod_NSB_calc_nano_100000", "step_mod_NSB_calc_nano_200000", "step_mod_NSB_calc_nano_500000", "step_mod_NSB_calc_nano_1000000")
all_foram <- c("step_mod_NSB_forams_100000", "step_mod_NSB_forams_200000", "step_mod_NSB_forams_500000", "step_mod_NSB_forams_1000000")
all_radio <- c("step_mod_NSB_radiolarians_100000", "step_mod_NSB_radiolarians_200000", "step_mod_NSB_radiolarians_500000", "step_mod_NSB_radiolarians_1000000")
all_diatom <- c("step_mod_NSB_diatoms_100000", "step_mod_NSB_diatoms_200000", "step_mod_NSB_diatoms_500000", "step_mod_NSB_diatoms_1000000")

###############################################

# And other info

org_list<-c("all_foram", "all_calc_nano", "all_radio", "all_diatom")

names <- c( "foraminifera", "calcareous_nanofossils", "radiolarians", "diatoms")

sizes <- c("100KY", "200KY", "500KY", "1000KY")

###############################################

# Calculating ratios of estimator coefficient (delta1:occupancy)

ratios <- data.frame(matrix(nrow = 16, ncol = 6))

colnames(ratios) <- c("data_set", "group", "binsize", "occupancy_estimator", "delta_estimator", "ratio_occupancy_delta")

count <- 0

for (j in 1:length(org_list)){
  
  org<-org_list[j]
  
  list <- get(org)
  
  for (i in 1:length(list)){
    
    count <- count +1
    
    name<- list[i]
    
    tmp<-get(list[i])
    
    occ<-as.numeric(tmp$coefficients[2])
    delta<-as.numeric(tmp$coefficients[3])
    ratio <- delta/occ
    
    ratios[count, 1]<-name
    ratios[count, 2]<-names[j]
    ratios[count, 3]<-sizes[i]
    ratios[count, 4]<-occ
    ratios[count, 5]<-delta
    ratios[count, 6]<-ratio
    
    
  }
  
}

foram_ratios <- ratios$ratio_occupancy_delta[ratios$group=="foraminifera"]
radio_ratios <- ratios$ratio_occupancy_delta[ratios$group=="radiolarians"]
calc_nano_ratios <- ratios$ratio_occupancy_delta[ratios$group=="calcareous_nanofossils"]
diatom_ratios <- ratios$ratio_occupancy_delta[ratios$group=="diatoms"]
