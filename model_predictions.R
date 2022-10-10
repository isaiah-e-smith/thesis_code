# fitting models and making predictions on extant organisms

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

############################3
# first, NSB data

files <- c("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_foraminifera_guest_2022-08-23_17-00-16.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_calc_nano_guest_2022-08-23_17-03-20.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_radiolarians_guest_2022-08-23_16-55-46.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_diatoms_guest_2022-08-23_17-04-54.tsv")

group_names <- c("forams", "calc_nano", "radiolarians", "diatoms")

bin_sizes <- c(1000000)

bin_sizes_strings<-c("1000000")

# create icosahedral grid, larger number = smaller grid cells
gr <- hexagrid(20, sp=TRUE)

# create empty data frame to fill with percent deviance values 
deviance_explained_delta1 <- data.frame(matrix(ncol = 3, nrow = (length(files)*length(bin_sizes)*2)))
colnames(deviance_explained_delta1) <- c("group", "bin_size", "value")



#########################################################################
# 1. Get the data
# load in data
#dat <- read.table("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/NSB_foraminifera_all_24_2_2022.tsv", sep = "\t", row.names=NULL, header=T)

dev_df_count <- 0


# select organism here
fi <- 1

ori_dat <- read.table(paste(files[fi]), sep = "\t", row.names=NULL, header=T, fill = T)


max(ori_dat$Age..Ma._Ogg.2020, na.rm=T)
length(ori_dat$gen_sp)
length(unique(ori_dat$gen_sp))

#########################################################################
# 2. taxonomic filtering


# Let's make sure every occurrence has a genus now
sum(is.na(ori_dat$Genus))
# Looking good (no NA's)

# Omit non-informative Genus entries
ori_dat <- ori_dat[ori_dat$Genus!="", ]
ori_dat <- ori_dat[is.na(ori_dat$Genus)==F, ]

# Let's make sure every occurrence has a species 
sum(is.na(ori_dat$Species))
# Looking good (no NA's)

# Omit non-informative species entries
ori_dat <- ori_dat[ori_dat$Species!="", ]
ori_dat <- ori_dat[is.na(ori_dat$Species)==F, ]

# Omit non-informative paleo-coordinate entries
ori_dat <- ori_dat[ori_dat$Estimated.Paleo.Latitude!="", ]
ori_dat <- ori_dat[ori_dat$Estimated.Paleo.Longitude!="", ]
ori_dat <- ori_dat[is.na(ori_dat$Estimated.Paleo.Latitude)==F, ]
ori_dat <- ori_dat[is.na(ori_dat$Estimated.Paleo.Longitude)==F, ]




# For now, keep everything


#########################################

# I will make a genus_species column, in case a species name is shared between two genera:

ori_dat$gen_sp <- paste(ori_dat$Genus, ori_dat$Species, sep = "_")

head(ori_dat$gen_sp)

length(ori_dat$gen_sp)

length(unique(ori_dat$gen_sp))

#set this column to whatever we want to do analysis on (genus or species), so we don't have to change the rest of the code
ori_dat$taxon <- ori_dat$gen_sp

head(ori_dat$taxon)



# Exclude occurrences with unassigned ages 
ori_dat <- subset(ori_dat, is.na(Age..Ma..Gradstein.et.al..2012)==F)

length(ori_dat$gen_sp)

length(unique(ori_dat$gen_sp))


# Taxonomic subset 
# Training: already extinct
# Testing: appeared in the last million years

extant<-c()
extinct<-c()

for (i in 1:length(unique(ori_dat$taxon))){
  un_tax<-unique(ori_dat$taxon)[i] 
  un_tax_ages<-ori_dat$Age..Ma..Gradstein.et.al..2012[ori_dat$taxon == un_tax]
  min_age<-min(un_tax_ages, na.rm = T)
  if(min_age <= (bin_sizes/1000000)){
    extant <- append(extant, un_tax)
  }else{
    extinct <- append(extinct, un_tax)
  }
}


# split into two data sets

train_ori_dat <- ori_dat[ori_dat$taxon %in% extinct,]
test_ori_dat <- ori_dat[ori_dat$taxon %in% extant,]

# is everything accounted for?
nrow(ori_dat) == nrow(train_ori_dat) + nrow(test_ori_dat)


# first, fit model with training set:
# Assign each occurrence to a bin

for (size in 1:length(bin_sizes)){
  
  dev_df_count <- dev_df_count + 1
  
  # For now, I am using 100,000 Year bins
  binsize <- bin_sizes[size]
  
  nb=ceiling(max(train_ori_dat$Age..Ma..Gradstein.et.al..2012)/(binsize/1000000))
  # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
  
  
  labs <- c(1:nb)
  
  labs <- as.character(labs)
  
  train_ori_dat$bin <- as.numeric(cut(train_ori_dat$Age..Ma..Gradstein.et.al..2012, breaks=(nb), labels = labs))
  
  
  # As of now, the more recent bins have smaller numbers; 
  # i.e. as we progress forward through time, the bins get smaller (the present is bin 1)
  
  
  # Identify two-timers
  # Which taxa are recorded in at least two bins
  
  t <- unique(subset(train_ori_dat, select=c(taxon, bin)))
  
  t.2 <- tapply(t$bin, t$taxon, length)
  
  length(t.2)
  
  t.2 <- t.2[t.2 > 1]
  
  length(t.2)
  
  t2.names <- names(t.2)
  
  t2.mat <- data.frame(subset(train_ori_dat, taxon %in% t2.names)) # analysis matrix (dataframe)
  
  # Remove removed factors
  #t2.mat$taxon <- factor(t2.mat$taxon)
  
  t2.inv <- unique(subset(t2.mat, select=c(taxon, bin)))
  
  
  # stratigraphic ranges of taxa
  orig <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, max)
  ext <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, min)
  
  
  # locate paleolatitude and paleolongitude coordinates on grid
  t2.mat$cells <- locate(gr, t2.mat[, c("Estimated.Paleo.Longitude", "Estimated.Paleo.Latitude")])
  
  
  ############################################################
  # this is where I begin modifying the code:
  
  
  
  # make empty data frame
  un_factor <- data.frame(matrix(ncol = 6, nrow = nrow(t2.mat)))
  
  
  
  # For calculating raw in-bin proportion (without considering cells sampled in previous bin)
  # and with "bridging gap" for missing intervals
  
  count <- 0
  
  # for delta1 versus delta2, just change binlist subset from "j +1" to "j + 2" and vice versa (two places to change)
  
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
      prev_bin <- binlist[j+1]
      
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
      
      # find raw proportion for this taxon in previously occupied bin
      if (is.na(un_before)==T | is.na(b)==T | length(un_before)==0 | length(b)==0){
        rpb <- NA
      }else{
        rpb <- length(un_before)/length(b)
      }
      
      
      
      prop <- (rp)/(rpb)
      
      if (is.infinite(prop)==T) {
        ans <- NA
      } else if (is.na(prop) == T){
        ans <- NA
      } else if (prop==0){
        ans <- NA 
        # or maybe leave 0? Shouldnt matter, since log(0) is still NA 
        # but a zero in the numerator would still represent decrease
        # maybe temporarily replace with 0.99 (smallest possible decrease)
        # shouldn't matter here, because there should not be a zero in the numerator
        # since we are skipping empty bins for each taxon
        
      } else {
        ans <- prop
      }
      
      
      # this is the delta value, which can be saved in the row for each taxon/bin combo
      lans <- log(ans)
      
      # save to data frame
      un_factor[count, 1] <- bin
      
      un_factor[count, 2] <- tx
      
      un_factor[count, 3] <- lans
      
      un_factor[count, 4] <- rp
      
      un_factor[count, 6] <- length(un)
      
    }
    
  }
  
  
  # name columns
  colnames(un_factor) <- c("bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled")
  
  # trim to only include rows that were filled
  un_factor <- un_factor[1:count,]
  
  
  
  # Fill Extinction column (1 = goes extinct or 0 = does not go extinct)
  
  # start by filling the Extinction column with all 0's (default: not going extinct in given bin)
  for (i in 1:nrow(un_factor)) {
    un_factor$extinction <- 0
  }
  
  # Now, for last occurrences (extinctions), assign Extinction value of 1
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
  
  un_factor_trim<-un_factor_trim[is.na(un_factor_trim$delta_1) == F,]
  
  # Save processed data 
  assign(paste("un_factor_trim", "NSB", group_names[fi], bin_sizes_strings[size], sep = "_"), (un_factor_trim))
  
  
  model1 <- glm(extinction ~ raw_prop_sampled * delta_1, 
                family = binomial(link ="logit"), data = un_factor_trim)
  
  step_model1 <- step(model1)
  
  # Save model output
  assign(paste("step_mod", "NSB", group_names[fi], bin_sizes_strings[size], sep = "_"), (step_model1))
  


}








# then, process testing set:

# Assign each occurrence to a bin

for (size in 1:length(bin_sizes)){
  
  dev_df_count <- dev_df_count + 1
  
  # For now, I am using 100,000 Year bins
  binsize <- bin_sizes[size]
  
  nb=ceiling(max(test_ori_dat$Age..Ma..Gradstein.et.al..2012)/(binsize/1000000))
  # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
  
  
  labs <- c(1:nb)
  
  labs <- as.character(labs)
  
  test_ori_dat$bin <- as.numeric(cut(test_ori_dat$Age..Ma..Gradstein.et.al..2012, breaks=(nb), labels = labs))
  
  
  # As of now, the more recent bins have smaller numbers; 
  # i.e. as we progress forward through time, the bins get smaller (the present is bin 1)
  
  
  # Identify two-timers
  # Which taxa are recorded in at least two bins
  
  t <- unique(subset(test_ori_dat, select=c(taxon, bin)))
  
  t.2 <- tapply(t$bin, t$taxon, length)
  
  length(t.2)
  
  t.2 <- t.2[t.2 > 1]
  
  length(t.2)
  
  t2.names <- names(t.2)
  
  t2.mat <- data.frame(subset(test_ori_dat, taxon %in% t2.names)) # analysis matrix (dataframe)
  
  # Remove removed factors
  #t2.mat$taxon <- factor(t2.mat$taxon)
  
  t2.inv <- unique(subset(t2.mat, select=c(taxon, bin)))
  
  
  # stratigraphic ranges of taxa
  orig <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, max)
  ext <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, min)
  
  
  # locate paleolatitude and paleolongitude coordinates on grid
  t2.mat$cells <- locate(gr, t2.mat[, c("Estimated.Paleo.Longitude", "Estimated.Paleo.Latitude")])
  
  
  ############################################################
  # this is where I begin modifying the code:
  
  
  
  # make empty data frame
  un_factor <- data.frame(matrix(ncol = 6, nrow = nrow(t2.mat)))
  
  
  
  # For calculating raw in-bin proportion (without considering cells sampled in previous bin)
  # and with "bridging gap" for missing intervals
  
  count <- 0
  
  # for delta1 versus delta2, just change binlist subset from "j +1" to "j + 2" and vice versa (two places to change)
  
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
      prev_bin <- binlist[j+1]
      
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
      
      # find raw proportion for this taxon in previously occupied bin
      if (is.na(un_before)==T | is.na(b)==T | length(un_before)==0 | length(b)==0){
        rpb <- NA
      }else{
        rpb <- length(un_before)/length(b)
      }
      
      
      
      prop <- (rp)/(rpb)
      
      if (is.infinite(prop)==T) {
        ans <- NA
      } else if (is.na(prop) == T){
        ans <- NA
      } else if (prop==0){
        ans <- NA 
        # or maybe leave 0? Shouldnt matter, since log(0) is still NA 
        # but a zero in the numerator would still represent decrease
        # maybe temporarily replace with 0.99 (smallest possible decrease)
        # shouldn't matter here, because there should not be a zero in the numerator
        # since we are skipping empty bins for each taxon
        
      } else {
        ans <- prop
      }
      
      
      # this is the delta value, which can be saved in the row for each taxon/bin combo
      lans <- log(ans)
      
      # save to data frame
      un_factor[count, 1] <- bin
      
      un_factor[count, 2] <- tx
      
      un_factor[count, 3] <- lans
      
      un_factor[count, 4] <- rp
      
      un_factor[count, 6] <- length(un)
      
    }
    
  }
  
  
  # name columns
  colnames(un_factor) <- c("bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled")
  
  # trim to only include rows that were filled
  un_factor <- un_factor[1:count,]
  
  
  
  # Fill Extinction column (1 = goes extinct or 0 = does not go extinct)
  
  # start by filling the Extinction column with all 0's (default: not going extinct in given bin)
  # for the test set, they all should still be extant (so 0)
  
  un_factor$extinction <- 0 
  
  head(un_factor)

}
  
#########################################
# save the processed data to unique variable
  
  
un_factor_test <- un_factor
  

  
# only look at the most recent bin (present bin) for extant organisms
un_factor_test_modern <-   un_factor_test[un_factor_test$bin == 1,]


# and make predictions (only one of the following chunks of code will work, based on which data set you just ran)

pred<-predict(step_mod_NSB_forams_1000000, newdata = un_factor_test_modern, type="response")
mean(pred, na.rm=T)
sd(pred, na.rm=T)
median(pred, na.rm=T)

par(mfrow = c(1,1), mar = (c(5,5,5,5)))
hist(pred, main = "NSB Foraminifera", xlab = "P(Ex)",
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7,
     col = "#487192")

pred<-predict(step_mod_NSB_calc_nano_1000000, newdata = un_factor_test_modern, type="response")
mean(pred, na.rm=T)
sd(pred, na.rm=T)
median(pred, na.rm=T)

par(mfrow = c(1,1), mar = (c(5,5,5,5)))
hist(pred, main = "NSB Calcareous Nanofossils", xlab = "P(Ex)",
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7,
     col = "#487192")

pred<-predict(step_mod_NSB_radiolarians_1000000, newdata = un_factor_test_modern, type="response")
mean(pred, na.rm=T)
sd(pred, na.rm=T)
median(pred, na.rm=T)

par(mfrow = c(1,1), mar = (c(5,5,5,5)))
hist(pred, main = "NSB Radiolarians", xlab = "P(Ex)",
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7,
     col = "#487192")

pred<-predict(step_mod_NSB_diatoms_1000000, newdata = un_factor_test_modern, type="response")
mean(pred, na.rm=T)
sd(pred, na.rm=T)
median(pred, na.rm=T)

par(mfrow = c(1,1), mar = (c(5,5,5,5)))
hist(pred, main = "NSB Diatoms", xlab = "P(Ex)",
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7,
     col = "#487192")





##########################3
# next, Triton data set

rm(list=ls())

load("~/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/2022-08-23_triton_RData.rdata")


group_names <- c("forams")

bin_sizes <- c(1000000)

bin_sizes_strings<-c("1000000")

# create icosahedral grid, larger number = smaller grid cells
gr <- hexagrid(20, sp=TRUE)

# create empty data frame to fill with percent deviance values 
deviance_explained_delta1 <- data.frame(matrix(ncol = 3, nrow = (length(files)*length(bin_sizes)*2)))
colnames(deviance_explained_delta1) <- c("group", "bin_size", "value")



#########################################################################
# 1. Get the data
# load in data
#dat <- read.table("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/NSB_foraminifera_all_24_2_2022.tsv", sep = "\t", row.names=NULL, header=T)

dev_df_count <- 0



ori_dat <- triton.pres


max(ori_dat$age, na.rm=T)
length(ori_dat$species)
length(unique(ori_dat$species))

#########################################################################
# 2. taxonomic filtering


# Let's make sure every occurrence has a species
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




# For now, keep everything


#########################################


ori_dat$taxon <- ori_dat$species

head(ori_dat$taxon)



# Exclude occurrences with unassigned ages 
ori_dat <- subset(ori_dat, is.na(age)==F)

length(ori_dat$species)

length(unique(ori_dat$species))


# Taxonomic subset 
# Training: already extinct
# Testing: appeared in the last million years

extant<-c()
extinct<-c()

for (i in 1:length(unique(ori_dat$taxon))){
  un_tax<-unique(ori_dat$taxon)[i] 
  un_tax_ages<-ori_dat$age[ori_dat$taxon == un_tax]
  min_age<-min(un_tax_ages, na.rm = T)
  if(min_age <= (bin_sizes/1000000)){
    extant <- append(extant, un_tax)
  }else{
    extinct <- append(extinct, un_tax)
  }
}


# split into two data sets

train_ori_dat <- ori_dat[ori_dat$taxon %in% extinct,]
test_ori_dat <- ori_dat[ori_dat$taxon %in% extant,]

# is everything accounted for?
nrow(ori_dat) == nrow(train_ori_dat) + nrow(test_ori_dat)


# first, fit model with training set:
# Assign each occurrence to a bin

for (size in 1:length(bin_sizes)){
  
  dev_df_count <- dev_df_count + 1
  
  # For now, I am using 100,000 Year bins
  binsize <- bin_sizes[size]
  
  nb=ceiling(max(train_ori_dat$age)/(binsize/1000000))
  # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
  
  
  labs <- c(1:nb)
  
  labs <- as.character(labs)
  
  train_ori_dat$bin <- as.numeric(cut(train_ori_dat$age, breaks=(nb), labels = labs))
  
  
  # As of now, the more recent bins have smaller numbers; 
  # i.e. as we progress forward through time, the bins get smaller (the present is bin 1)
  
  
  # Identify two-timers
  # Which taxa are recorded in at least two bins
  
  t <- unique(subset(train_ori_dat, select=c(taxon, bin)))
  
  t.2 <- tapply(t$bin, t$taxon, length)
  
  length(t.2)
  
  t.2 <- t.2[t.2 > 1]
  
  length(t.2)
  
  t2.names <- names(t.2)
  
  t2.mat <- data.frame(subset(train_ori_dat, taxon %in% t2.names)) # analysis matrix (dataframe)
  
  # Remove removed factors
  #t2.mat$taxon <- factor(t2.mat$taxon)
  
  t2.inv <- unique(subset(t2.mat, select=c(taxon, bin)))
  
  
  # stratigraphic ranges of taxa
  orig <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, max)
  ext <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, min)
  
  
  # locate paleolatitude and paleolongitude coordinates on grid
  t2.mat$cells <- locate(gr, t2.mat[, c("pal.long", "pal.lat")])
  
  
  ############################################################
  # this is where I begin modifying the code:
  
  
  
  # make empty data frame
  un_factor <- data.frame(matrix(ncol = 6, nrow = nrow(t2.mat)))
  
  
  
  # For calculating raw in-bin proportion (without considering cells sampled in previous bin)
  # and with "bridging gap" for missing intervals
  
  count <- 0
  
  # for delta1 versus delta2, just change binlist subset from "j +1" to "j + 2" and vice versa (two places to change)
  
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
      prev_bin <- binlist[j+1]
      
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
      
      # find raw proportion for this taxon in previously occupied bin
      if (is.na(un_before)==T | is.na(b)==T | length(un_before)==0 | length(b)==0){
        rpb <- NA
      }else{
        rpb <- length(un_before)/length(b)
      }
      
      
      
      prop <- (rp)/(rpb)
      
      if (is.infinite(prop)==T) {
        ans <- NA
      } else if (is.na(prop) == T){
        ans <- NA
      } else if (prop==0){
        ans <- NA 
        # or maybe leave 0? Shouldnt matter, since log(0) is still NA 
        # but a zero in the numerator would still represent decrease
        # maybe temporarily replace with 0.99 (smallest possible decrease)
        # shouldn't matter here, because there should not be a zero in the numerator
        # since we are skipping empty bins for each taxon
        
      } else {
        ans <- prop
      }
      
      
      # this is the delta value, which can be saved in the row for each taxon/bin combo
      lans <- log(ans)
      
      # save to data frame
      un_factor[count, 1] <- bin
      
      un_factor[count, 2] <- tx
      
      un_factor[count, 3] <- lans
      
      un_factor[count, 4] <- rp
      
      un_factor[count, 6] <- length(un)
      
    }
    
  }
  
  
  # name columns
  colnames(un_factor) <- c("bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled")
  
  # trim to only include rows that were filled
  un_factor <- un_factor[1:count,]
  
  
  
  # Fill Extinction column (1 = goes extinct or 0 = does not go extinct)
  
  # start by filling the Extinction column with all 0's (default: not going extinct in given bin)
  for (i in 1:nrow(un_factor)) {
    un_factor$extinction <- 0
  }
  
  # Now, for last occurrences (extinctions), assign Extinction value of 1
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
  
  un_factor_trim<-un_factor_trim[is.na(un_factor_trim$delta_1) == F,]
  
  
  model1 <- glm(extinction ~ raw_prop_sampled * delta_1, 
                family = binomial(link ="logit"), data = un_factor_trim)
  
  step_model1 <- step(model1)
  
  # Save model output
  assign(paste("step_mod", "triton", bin_sizes_strings[size], sep = "_"), (step_model1))
  

  
}



# then, process testing set:

# Assign each occurrence to a bin

for (size in 1:length(bin_sizes)){
  
  dev_df_count <- dev_df_count + 1
  
  # For now, I am using 100,000 Year bins
  binsize <- bin_sizes[size]
  
  nb=ceiling(max(test_ori_dat$age)/(binsize/1000000))
  # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
  
  
  labs <- c(1:nb)
  
  labs <- as.character(labs)
  
  test_ori_dat$bin <- as.numeric(cut(test_ori_dat$age, breaks=(nb), labels = labs))
  
  
  # As of now, the more recent bins have smaller numbers; 
  # i.e. as we progress forward through time, the bins get smaller (the present is bin 1)
  
  
  # Identify two-timers
  # Which taxa are recorded in at least two bins
  
  t <- unique(subset(test_ori_dat, select=c(taxon, bin)))
  
  t.2 <- tapply(t$bin, t$taxon, length)
  
  length(t.2)
  
  t.2 <- t.2[t.2 > 1]
  
  length(t.2)
  
  t2.names <- names(t.2)
  
  t2.mat <- data.frame(subset(test_ori_dat, taxon %in% t2.names)) # analysis matrix (dataframe)
  
  # Remove removed factors
  #t2.mat$taxon <- factor(t2.mat$taxon)
  
  t2.inv <- unique(subset(t2.mat, select=c(taxon, bin)))
  
  
  # stratigraphic ranges of taxa
  orig <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, max)
  ext <-  tapply(as.numeric(t2.inv$bin), t2.inv$taxon, min)
  
  
  # locate paleolatitude and paleolongitude coordinates on grid
  t2.mat$cells <- locate(gr, t2.mat[, c("pal.long", "pal.lat")])
  
  
  ############################################################
  # this is where I begin modifying the code:
  
  
  
  # make empty data frame
  un_factor <- data.frame(matrix(ncol = 6, nrow = nrow(t2.mat)))
  
  
  
  # For calculating raw in-bin proportion (without considering cells sampled in previous bin)
  # and with "bridging gap" for missing intervals
  
  count <- 0
  
  # for delta1 versus delta2, just change binlist subset from "j +1" to "j + 2" and vice versa (two places to change)
  
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
      prev_bin <- binlist[j+1]
      
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
      
      # find raw proportion for this taxon in previously occupied bin
      if (is.na(un_before)==T | is.na(b)==T | length(un_before)==0 | length(b)==0){
        rpb <- NA
      }else{
        rpb <- length(un_before)/length(b)
      }
      
      
      
      prop <- (rp)/(rpb)
      
      if (is.infinite(prop)==T) {
        ans <- NA
      } else if (is.na(prop) == T){
        ans <- NA
      } else if (prop==0){
        ans <- NA 
        # or maybe leave 0? Shouldnt matter, since log(0) is still NA 
        # but a zero in the numerator would still represent decrease
        # maybe temporarily replace with 0.99 (smallest possible decrease)
        # shouldn't matter here, because there should not be a zero in the numerator
        # since we are skipping empty bins for each taxon
        
      } else {
        ans <- prop
      }
      
      
      # this is the delta value, which can be saved in the row for each taxon/bin combo
      lans <- log(ans)
      
      # save to data frame
      un_factor[count, 1] <- bin
      
      un_factor[count, 2] <- tx
      
      un_factor[count, 3] <- lans
      
      un_factor[count, 4] <- rp
      
      un_factor[count, 6] <- length(un)
      
    }
    
  }
  
  
  # name columns
  colnames(un_factor) <- c("bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled")
  
  # trim to only include rows that were filled
  un_factor <- un_factor[1:count,]
  
  
  
  # Fill Extinction column (1 = goes extinct or 0 = does not go extinct)
  
  # start by filling the Extinction column with all 0's (default: not going extinct in given bin)
  # for the test set, they all should still be extant (so 0)
  
  un_factor$extinction <- 0 
  
  head(un_factor)
  
}

#########################################
# save the processed data to unique variable

un_factor_test <- un_factor


# only keep most recent bin (present bin)
un_factor_test_modern <-   un_factor_test[un_factor_test$bin == 1,]



#### whole data set averages of probabilities

pred<-predict(step_mod_triton_1000000, newdata = un_factor_test_modern, type="response")
mean(pred, na.rm=T)
sd(pred, na.rm=T)
median(pred, na.rm=T)

par(mfrow = c(1,1), mar = (c(5,5,5,5)))
hist(pred, main = "Triton", xlab = "P(Ex)",
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7,
     col = "#F0BED5")













# the following code only works for the NSB foraminifera data set,
# and was used for examples of specific organisms

# put it all together (for NSB foraminifera):

par( mfrow= c(2,3), mai = c(.6,.6,.6,.2) )

## 1
this_tax <- "Globigerinoides_extremus"
newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(13),
     y = 0.4,
     labels = paste('P(Ex):', sprintf('%.2f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)

## 2


this_tax <- "Pulleniatina_primalis"
newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(12),
     y = 0.305,
     labels = paste('P(Ex):', sprintf('%.2f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)




#3

this_tax <- "Subbotina_corpulenta"

newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     ylim = c(0,0.1),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(9),
     y = 0.079,
     labels = paste('P(Ex):', sprintf('%.2f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)







#4


this_tax <- "Hastigerina_siphonifera"
newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(12),
     y = 0.08,
     labels = paste('P(Ex):', sprintf('%.2f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)



#5

this_tax <- "Globorotalia_crassaformis"
newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(12),
     y = 0.54,
     labels = paste('P(Ex):', sprintf('%.4f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)





#6


this_tax <- "Globigerinoides_ruber"
newdata = un_factor_test_modern[un_factor_test_modern$taxon == this_tax,]
pred<-predict(step_mod_NSB_forams_1000000, newdata = newdata, type="response")

plot_dat <-  un_factor_test[un_factor_test$taxon == this_tax,]
plot(x = plot_dat$bin, 
     y = plot_dat$raw_prop_sampled, 
     type = "step",
     xlim = c(30, 0),
     col = c("darkred"),
     ylab = "Proportion Sampled",
     xlab = "MYA",
     main = sub("_"," ", this_tax),
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.7)

text(x = c(12),
     y = 0.5,
     labels = paste('P(Ex):', sprintf('%.4f',pred)),
     xpd = NA,
     srt = 0,
     adj = 0.95,
     cex = 1.8)



dev.off()
