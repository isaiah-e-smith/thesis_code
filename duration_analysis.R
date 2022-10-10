# Calculating completeness metrics of different data sets 

library(divDyn)

files <- c("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_foraminifera_guest_2022-08-23_17-00-16.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_calc_nano_guest_2022-08-23_17-03-20.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_radiolarians_guest_2022-08-23_16-55-46.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_diatoms_guest_2022-08-23_17-04-54.tsv")

group_names <- c("forams", "calc_nano", "radiolarians","diatoms")

bin_sizes <- c(100000, 200000, 500000, 1000000)

bin_sizes_strings<-c("100000", "200000", "500000", "1000000")


# prepare empty dataframe to be filled
duration_data <- data.frame(matrix( ncol = 3, nrow = 10000))

colnames(duration_data) <- c("group", "duration", "present_in_recent")


dev_df_count <- 0

for (fi in 1:length(files)){ 
  
  ori_dat <- read.table(paste(files[fi]), sep = "\t", row.names=NULL, header=T, fill = T)
  
  
  #########################################################################
  # 2. taxonomic filtering
  
  
  # Let's make sure every occurrence has a genus 
  sum(is.na(ori_dat$Genus))
  # Looking good (no NA's)
  
  # Omit non-informative genus entries
  ori_dat <- ori_dat[ori_dat$Genus!="", ]
  ori_dat <- ori_dat[is.na(ori_dat$Genus)==F, ]
  
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
  # I suppose this is redundant, since we drop the most recent bin from model analysis anyways
  ori_dat <- subset(ori_dat, is.na(Age..Ma..Gradstein.et.al..2012)==F)
  
  # And exclude occurrences with ages of 0 (still extant)
  ori_dat <- ori_dat[ori_dat$Age..Ma..Gradstein.et.al..2012!="0", ]
  
  #####
  
  for (n in 1:length(unique(ori_dat$gen_sp))){
    
    dev_df_count <- dev_df_count + 1
    
    ttax <- unique(ori_dat$gen_sp)[n]
    
    max_age <- max(ori_dat$Age..Ma..Gradstein.et.al..2012[ori_dat$gen_sp == ttax])
    
    min_age <- min(ori_dat$Age..Ma..Gradstein.et.al..2012[ori_dat$gen_sp == ttax])
    
    duration <- max_age - min_age
    
    duration_data$group[dev_df_count] <- (group_names[fi])
    
    duration_data$duration[dev_df_count] <- duration
    
    if (min_age <= 0.1){ #the value specified here can be changed to specify how recent can an occurrence be for a species to be considered "present in The Recent"
      duration_data$present_in_recent[dev_df_count] <- 1
    } else{
      duration_data$present_in_recent[dev_df_count] <- 0
    }
    
    
  }
  
} 

table(duration_data$present_in_recent[duration_data$group == "forams"])
table(duration_data$present_in_recent[duration_data$group == "calc_nano"])
table(duration_data$present_in_recent[duration_data$group == "radiolarians"])
table(duration_data$present_in_recent[duration_data$group == "diatoms"])


mean(duration_data$duration[duration_data$group == "forams"], na.rm = T)
mean(duration_data$duration[duration_data$group == "calc_nano"], na.rm = T)
mean(duration_data$duration[duration_data$group == "radiolarians"], na.rm = T)
mean(duration_data$duration[duration_data$group == "diatoms"], na.rm = T)

sd(duration_data$duration[duration_data$group == "forams"], na.rm = T)
sd(duration_data$duration[duration_data$group == "calc_nano"], na.rm = T)
sd(duration_data$duration[duration_data$group == "radiolarians"], na.rm = T)
sd(duration_data$duration[duration_data$group == "diatoms"], na.rm = T)

median(duration_data$duration[duration_data$group == "forams"], na.rm = T)
median(duration_data$duration[duration_data$group == "calc_nano"], na.rm = T)
median(duration_data$duration[duration_data$group == "radiolarians"], na.rm = T)
median(duration_data$duration[duration_data$group == "diatoms"], na.rm = T)

par( mfrow= c(2,3), mai = c(.6,.6,.2,0.1) )

hist(duration_data$duration[duration_data$group == "forams"], 
     xlab = "Duration of Species (Millions of Years)",
     main = "NSB Foraminifera",
     ylim = c(0,1000),
     col = "#487192",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.main = 1.5)

hist(duration_data$duration[duration_data$group == "calc_nano"], 
     xlab = "Duration of Species (Millions of Years)",
     main = "NSB Calcareous Nanofossils",
     ylim = c(0,1000),
     col = "#487192",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.main = 1.5)

hist(duration_data$duration[duration_data$group == "radiolarians"], 
     xlab = "Duration of Species (Millions of Years)",
     main = "NSB Radiolarians",
     ylim = c(0,1000),
     col = "#487192",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.main = 1.5)

hist(duration_data$duration[duration_data$group == "diatoms"], 
     xlab = "Duration of Species (Millions of Years)",
     main = "NSB Diatoms",
     ylim = c(0,1000),
     col = "#487192",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.main = 1.5)

