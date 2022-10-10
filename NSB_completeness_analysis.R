# Calculating completeness metrics of Neptune Sandbox Berlin data sets

library(divDyn)

files <- c("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_foraminifera_guest_2022-08-23_17-00-16.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_calc_nano_guest_2022-08-23_17-03-20.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_radiolarians_guest_2022-08-23_16-55-46.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_diatoms_guest_2022-08-23_17-04-54.tsv")

group_names <- c("forams", "calc_nano", "radiolarians","diatoms")

bin_sizes <- c(100000, 200000, 500000, 1000000)

bin_sizes_strings<-c("100000", "200000", "500000", "1000000")


# prepare empty dataframe to be filled
completeness_data <- data.frame(matrix(nrow = length(bin_sizes)*length(files), ncol = 6))

colnames(completeness_data) <- c("group", "binsize", "bin_avg_proportion_part_timer", "all_bin_avg_three_timer","three_timer_SD", "whole_dataset_SCM")


dev_df_count <- 0

for (fi in 1:length(files)){ 
  
  ori_dat <- read.table(paste(files[fi]), sep = "\t", row.names=NULL, header=T, fill = T)
  
  
  #########################################################################
  # 2. taxonomic filtering
  
  
  # Let's make sure every occurrence has a genus 
  sum(is.na(ori_dat$Genus))
  # Looking good (no NA's)
  
  # Omit non-informative Genus entries
  ori_dat <- ori_dat[ori_dat$Genus!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$Genus)==F, ]
  
  # Let's make sure every occurrence has a species 
  sum(is.na(ori_dat$Species))
  # Looking good (no NA's)
  
  # Omit non-informative Genus entries
  ori_dat <- ori_dat[ori_dat$Species!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$Species)==F, ]
  
  # Omit non-informative paleo-coordinate entries
  ori_dat <- ori_dat[ori_dat$Estimated.Paleo.Latitude!="", ]
  ori_dat <- ori_dat[ori_dat$Estimated.Paleo.Longitude!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$Estimated.Paleo.Latitude)==F, ]
  ori_dat <- ori_dat[is.na(ori_dat$Estimated.Paleo.Longitude)==F, ]
  
  
  # Taxonomic subset 
  
  # For now, keep everything
  
  
  #########################################
  
  # I will make a genus_species column, in case a species name is shared between two genera:
  
  ori_dat$gen_sp <- paste(ori_dat$Genus, ori_dat$Species, sep = "_")
  
  head(ori_dat$gen_sp)
  
  #set this column to whatever we want to do analysis on (genus or species), so we don't have to change the rest of the code
  ori_dat$taxon <- ori_dat$gen_sp
  #ori_dat$taxon <- ori_dat$Genus
  
  head(ori_dat$taxon)
  
  # Exclude occurrences with unassigned ages 
  ori_dat <- subset(ori_dat, is.na(Age..Ma..Gradstein.et.al..2012)==F)
  
  # Assign each occurrence to a bin
  
  for (size in 1:length(bin_sizes)){
    
    dev_df_count <- dev_df_count + 1
    
    print(dev_df_count)
    
    # select bin size
    binsize <- bin_sizes[size]
    
    nb=ceiling(max(ori_dat$Age..Ma..Gradstein.et.al..2012)/(binsize/1000000))
    # "(binsize/1000000)" accounts for the fact that time is reported in millions of years
    
    
    labs <- c(1:nb)
    
    labs <- as.character(labs)
    
    ori_dat$bin <- as.numeric(cut(ori_dat$Age..Ma..Gradstein.et.al..2012, breaks=(nb), labels = labs))
    
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


############################################## 
# Plotting SCM

par(mfrow = c(1,4), mar = c(3, 5, 5, 2)) #mar=c(0.5,5,5,0.1)

# first row, 1 binsize
barplot(c(completeness_data$whole_dataset_SCM[1], completeness_data$whole_dataset_SCM[5], completeness_data$whole_dataset_SCM[9], 
          completeness_data$whole_dataset_SCM[13]), ylim = c(0,1),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
        col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
        cex.names = 1.2, cex.lab = 1.5, las = 0)

title(ylab = "Binsize: 1e+5", line=3, cex.lab=1.5)


# second row, next binsize
barplot(c(completeness_data$whole_dataset_SCM[2], completeness_data$whole_dataset_SCM[6], completeness_data$whole_dataset_SCM[10], 
          completeness_data$whole_dataset_SCM[14]), ylim = c(0,1),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
        col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
        cex.names = 1.2, cex.lab = 1.5, las = 0)

title(ylab = "Binsize: 2e+5", line=3, cex.lab=1.5)


# third row, next binsize
barplot(c(completeness_data$whole_dataset_SCM[3], completeness_data$whole_dataset_SCM[7], completeness_data$whole_dataset_SCM[11], 
          completeness_data$whole_dataset_SCM[15]), ylim = c(0,1),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
        col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
        cex.names = 1.2, cex.lab = 1.5, las = 0)

title(ylab = "Binsize: 5e+5", line=3, cex.lab=1.5)


# fourth row, next binsize
barplot(c(completeness_data$whole_dataset_SCM[4], completeness_data$whole_dataset_SCM[8], completeness_data$whole_dataset_SCM[12], 
          completeness_data$whole_dataset_SCM[16]), ylim = c(0,1),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
        col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
        cex.names = 1.2, cex.lab = 1.5, las = 0)

title(ylab = "Binsize: 1e+6", line=3, cex.lab=1.5)


# Add title to whole multi-plot
mtext("Simple Completeness Metric", side = 3, line = -3, outer = TRUE, cex = 1.6)


############################################## 
# Plotting 3-timer completeness

par(mfrow = c(1,4), mar = c(3, 5, 5, 2)) #mar=c(0.5,5,5,0.1)

# first binsize

b1<- c(completeness_data$all_bin_avg_three_timer[1], completeness_data$all_bin_avg_three_timer[5], completeness_data$all_bin_avg_three_timer[9], 
        completeness_data$all_bin_avg_three_timer[13])

b1sd <- c(completeness_data$three_timer_SD[1], completeness_data$three_timer_SD[5], completeness_data$three_timer_SD[9], 
          completeness_data$three_timer_SD[13])

bp1<-barplot(b1, ylim = c(0,1.3),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
        col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
        cex.names = 1.2, cex.lab = 1.5, las = 0)

arrows(x0=bp1, y0=b1-b1sd, x1=bp1, y1=b1+b1sd,length = 0.05, angle = 90, code = 3)

title(ylab = "Binsize: 1e+5", line=3, cex.lab=1.5)


# second binsize

b1<- c(completeness_data$all_bin_avg_three_timer[2], completeness_data$all_bin_avg_three_timer[6], completeness_data$all_bin_avg_three_timer[10], 
       completeness_data$all_bin_avg_three_timer[14])

b1sd <- c(completeness_data$three_timer_SD[2], completeness_data$three_timer_SD[6], completeness_data$three_timer_SD[10], 
          completeness_data$three_timer_SD[14])

bp1<-barplot(b1, ylim = c(0,1.3),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
             col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
             cex.names = 1.2, cex.lab = 1.5, las = 0)

arrows(x0=bp1, y0=b1-b1sd, x1=bp1, y1=b1+b1sd,length = 0.05, angle = 90, code = 3)

title(ylab = "Binsize: 2e+5", line=3, cex.lab=1.5)



# third binsize

b1<- c(completeness_data$all_bin_avg_three_timer[3], completeness_data$all_bin_avg_three_timer[7], completeness_data$all_bin_avg_three_timer[11], 
       completeness_data$all_bin_avg_three_timer[15])

b1sd <- c(completeness_data$three_timer_SD[3], completeness_data$three_timer_SD[7], completeness_data$three_timer_SD[11], 
          completeness_data$three_timer_SD[15])

bp1<-barplot(b1, ylim = c(0,1.3),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
             col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
             cex.names = 1.2, cex.lab = 1.5, las = 0)

arrows(x0=bp1, y0=b1-b1sd, x1=bp1, y1=b1+b1sd,length = 0.05, angle = 90, code = 3)

title(ylab = "Binsize: 5e+5", line=3, cex.lab=1.5)



# fourth binsize

b1<- c(completeness_data$all_bin_avg_three_timer[4], completeness_data$all_bin_avg_three_timer[8], completeness_data$all_bin_avg_three_timer[12], 
       completeness_data$all_bin_avg_three_timer[16])

b1sd <- c(completeness_data$three_timer_SD[4], completeness_data$three_timer_SD[8], completeness_data$three_timer_SD[12], 
          completeness_data$three_timer_SD[16])

bp1<-barplot(b1, ylim = c(0,1.3),ylab = "", xlab = " ", cex.lab = 1.5, cex.axis = 1.3, 
             col = c("gray", "gray", "dark blue", "dark blue"), names.arg = c("F", "CN", "R", "D"), 
             cex.names = 1.2, cex.lab = 1.5, las = 0)

arrows(x0=bp1, y0=b1-b1sd, x1=bp1, y1=b1+b1sd,length = 0.05, angle = 90, code = 3)

title(ylab = "Binsize: 1e+6", line=3, cex.lab=1.5)

# Add title to whole multi-plot
mtext("Three-Timer Completeness, Averaged Across Bins", side = 3, line = -3, outer = TRUE, cex = 1.6)

# 






################################

mean(completeness_data$bin_avg_three_timer[completeness_data$group == "forams"])
sd(completeness_data$bin_avg_three_timer[completeness_data$group == "forams"])

mean(completeness_data$bin_avg_three_timer[completeness_data$group == "calc_nano"])
sd(completeness_data$bin_avg_three_timer[completeness_data$group == "calc_nano"])

mean(completeness_data$bin_avg_three_timer[completeness_data$group == "radiolarians"])
sd(completeness_data$bin_avg_three_timer[completeness_data$group == "radiolarians"])

mean(completeness_data$bin_avg_three_timer[completeness_data$group == "diatoms"])
sd(completeness_data$bin_avg_three_timer[completeness_data$group == "diatoms"])


mean(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "forams"])
sd(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "forams"])

mean(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "calc_nano"])
sd(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "calc_nano"])

mean(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "radiolarians"])
sd(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "radiolarians"])

mean(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "diatoms"])
sd(completeness_data$bin_avg_proportion_part_timer[completeness_data$group == "diatoms"])


