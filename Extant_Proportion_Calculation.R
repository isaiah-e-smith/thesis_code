# Proportion Extant


getwd()
setwd("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/")


files <- c("/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_foraminifera_guest_2022-08-23_17-00-16.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_calc_nano_guest_2022-08-23_17-03-20.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_radiolarians_guest_2022-08-23_16-55-46.tsv",
           "/Users/isaiahsmith/Desktop/Research_Projects/FAU_Thesis/Foraminifera_Occupancy_Trajectory/Thesis_Data_23_August_2022/NSB_diatoms_guest_2022-08-23_17-04-54.tsv")




#########################################################################
#NSB

fi <- 4 #select which data set

bin_size <- 1000000 #choose bin size

ori_dat <- read.table(paste(files[fi]), sep = "\t", row.names=NULL, header=T, fill = T)


# Exclude occurrences with unassigned ages 
ori_dat <- subset(ori_dat, is.na(Age..Ma..Gradstein.et.al..2012)==F)

ori_dat$gen_sp <- paste(ori_dat$Genus, ori_dat$Species, sep = "_")

ori_dat$taxon <- ori_dat$gen_sp

extant<-c()

extinct<-c()

for (i in 1:length(unique(ori_dat$taxon))){
  un_tax<-unique(ori_dat$taxon)[i] 
  un_tax_ages<-ori_dat$Age..Ma..Gradstein.et.al..2012[ori_dat$taxon == un_tax]
  min_age<-min(un_tax_ages, na.rm = T)
  if(min_age <= (bin_size/1000000)){
    extant <- append(extant, un_tax)
  }else{
    extinct <- append(extinct, un_tax)
  }
}

length(extant)/(length(extant)+length(extinct))


###########
# Triton

bin_size <- 1000000

dat <- triton.pres #object from downloadalbe Triton file

ori_dat <- dat

# Exclude occurrences with unassigned ages 
ori_dat <- subset(ori_dat, is.na(age)==F)

ori_dat$taxon <- ori_dat$species

extant<-c()

extinct<-c()

for (i in 1:length(unique(ori_dat$taxon))){
  un_tax<-unique(ori_dat$taxon)[i] 
  un_tax_ages<-ori_dat$age[ori_dat$taxon == un_tax]
  min_age<-min(un_tax_ages, na.rm = T)
  if(min_age <= (bin_size/1000000)){
    extant <- append(extant, un_tax)
  }else{
    extinct <- append(extinct, un_tax)
  }
}

length(extant)/(length(extant)+length(extinct))








