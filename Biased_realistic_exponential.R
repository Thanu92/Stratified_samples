
#Function to combine two rows in a dataframe
combine_rows <- function(data, row1, row2) {
  data[row2, ] <- data[row1, ] + data[row2, ]
  data[-row1, ]
}

#Read the fish dataset in R
fish <- read.table("final_alignment.phylip")

# install.packages("tidyverse")
library(tidyverse)
library(dplyr)
#install.packages("rfishbase")                   
library(rfishbase)
library(stats)
library(grid)
library(futile.logger)
library(VennDiagram)
library(dplyr)
library(MASS)

#Rename columns
colnames(fish) <- c("species_name","seq")

#Extract CO1 nucleotides (This step is to extract sequences with CO1 gene. So that we can remove sequences without CO1 gene)
fishCO1 <- substr(fish$seq,2292,2973)

#Check the class of fishCO1 dataset
class(fishCO1)

#Convert to dataframe
fishCO1 <- data.frame(fishCO1)

#Check the converted dataset
class(fishCO1)

#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(fish$species_name, fishCO1)

#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)

#Convert to dataframe
df_0 <- as.data.frame(df_0)

#Bind the column "species name" of fish dataset with new dataframe
df_1 <- cbind(fish$species_name,df_0)

#Remove rows with empty CO1 regions (For the phylogenetic placement we do need CO1 region. Therefore the empty CO1 regions are useless for further analysis)
df_2 <- subset(df_1, df_1$df_0 != "")

#Rename column names
colnames(df_2) <- c("species_name","seq")

#Keep only rows with sequences
fish_multigene <- merge(fish,df_2, by= "species_name")

#New dataframe with first two columns (Get the dataframe which we need for the down stream analysis)
fish_multigene <- fish_multigene[, 1:2]

#Rename columns
colnames(fish_multigene) <- c("species_name","seq")
class(fish_multigene)

#Read the fish taxonomy spread sheet in R
fish_taxonomy <- read.csv("/home/thanu/Desktop/FishData/Stratified/PFCtaxonomy .csv")

#Extract family column and genus_species column for downstream analysis
df_fish_family_species <- fish_taxonomy[,c("family","genus.species")]

#rename columns
colnames(df_fish_family_species) <- c("family","species_name")

#Replace the space with "_" and covert the resultant list into a df
df_fish_family_species1 <- as.data.frame(gsub(' ','_',df_fish_family_species$species_name))

#Bind the colum "family" of fish_family_species dataset with new dataframe df_fish_family_species1
newdf_fish_taxon <- cbind(df_fish_family_species1,df_fish_family_species$family)

#rename columns
colnames(newdf_fish_taxon) <- c("species_name","family")

class(newdf_fish_taxon)
#merge fish taxon and sequences
#df_intersect <- intersect(newdf_fish_taxon,fish_multigene)
df_fish <- merge(newdf_fish_taxon,fish_multigene,by="species_name")

#Checking the families with the number of species in each family and convert to a dataframe
df_tt <- data.frame(table(df_fish$family))

#Rename columns
colnames(df_tt) <- c("family","species_count")

#Check the umber of families in fish sample
length(unique(df_fish$family))

# Fishbase dataset matched against the information obtained from FishTreeofLife dataset 
# Extract all of the species names that are available on FishBase
AllFish <- fishbase

#Paste Genus and Species columns together to get the species name 
AllFish$FullName <- paste(AllFish$Genus, AllFish$Species) 
FishBaseSpecies <- AllFish$FullName # 33104 species names
#Match the species labels from FishTreeofLife dataset with the species names from FishBase

#Make FishBaseSpecies into a dataframe first so it can be merged with the family column of AllFish dataset

dfFishBaseSpecies <- data.frame(FishBaseSpecies)

#Checking the columns of fishBase dataset
names(fishbase)

#combine the species name and family columns together in fishbase data
FB_family_species <- cbind(AllFish$Family,dfFishBaseSpecies)

#Rename columns
colnames(FB_family_species) <- c("family","species_name")
class(FB_family_species)
#Checking the families with the number of species in each family
tt_FB <- table(FB_family_species$family)
df_tt_FB <- data.frame(tt_FB <- table(FB_family_species$family))

#Rename columns
colnames(df_tt_FB) <- c("family","species_count")

#Check the number of families in fish sample
length(unique(FB_family_species$family))

#Merge the two dataframes of rfish base and fishTreeofLife dataset by "family" column
df_merged <- merge(df_tt, df_tt_FB, by="family")

#Rename columns
colnames(df_merged) <- c("family","FishTreeofLife_species_count","FishBase_species_count")
#There are 364 familes common in fishTree and FishBase
#By looking at this dataset, I found Psychrolutidae has more species count in fishTree than rFishbase. 
#When checking the literature deeply, I was able to identify that Psychrolutidae is a sub set of Cottilae in old work.
#Likewise Alepisaundridae has more species count in fishTree than fishBase. 
#According to literature, there is a close relationship between Alepisauridae and Paralepididae.
#Usually, fishTree dataset should be a subset of Fishbase.
#Hence, this count shows (either the fishbase is incomplete or a descripency in taxonomy between 2 datasets)
#As the mentioned families are correlated, I have combined the two families together for down stream analysis

#this is to get meaningful values (which implies fishTree is a subset of FishBase)
#Combine Psychrolutidae and cottidae family group (row numbers for Psychrolutidae and cottidae are 249 and 102 respectively)
df_merged_com1 <- combine_rows(df_merged, 294, 102)

#Combine Alepisauridae and Paralepididae familes (row numbers for Psychrolutidae and cottidae are 10 and 257 respectively)
df_merged_com2 <- combine_rows(df_merged_com1, 10, 257)
tail(df_merged_com2)

#Calculation
#Get the species count percentage per family by dividing the FishTreeofLife species count by FishBase species count
df_merged_1<- df_merged_com2 %>% mutate(present_species_count_percentage = (FishTreeofLife_species_count/ FishBase_species_count)*100)

#df_merged_1
#missing data
df_merged_2<- df_merged_1 %>% mutate(missing_percentage=100-(FishTreeofLife_species_count/ FishBase_species_count)*100)

#Let's check a venn diagram,
# Chart
venn.diagram(
  x = list(df_fish$family,FB_family_species$family),
  fill = c("green", "red"),
  category.names = c("FishTree", "FishBase"),
  filename = 'venn_diagram',
  output=TRUE
)
#According to this 364 families are common in both datasets.
#14 familes are there only in fishTree and not in fishbase(Which cannot be acceptable according to our assumption: fishTree is a subset of FishBase).
#I checked these families deeply and found descripency of taxonomy between two datasets
#Hence, we can exclude these 14 for downstream analysis.
#185 familes can be seen in fishbase but not in fishTree (Which can be acceptable)
#So, I'm going to include these 185 families for downstream analysis.
#Now there are 547 families together for the missing data analysis 
#(present species percentage for these 185 families are 0% therefore, the missing species percentage is 100%)
#Number of families in the genetic dataset is 378 (in fish tree of life dataset).

#Union of the two dataset
unio <- union(df_fish$family,FB_family_species$family)
class(unio)

#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))

#Add the column the missing species percentage which is 100%
setdif_FishBase$missing_percentage <- 100

#563-378 = 185 (fish families in fishTree but not in fish base)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","missing_percentage")

#make a new dataframe by extracting family and missing percentage columns
df_merged_3 <- df_merged_2[,c("family","missing_percentage")]

#Adding 185 families with 100% missing species percentage to 362 families, 
#the percentages are zeros here for present species
new <- rbind(df_merged_3,setdif_FishBase)
class(new)
#checking the column names of the new dataframe
names(new)

#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))

#Add the column the present species percentage which is 0%
setdif_FishBase$present_species_count_percentage <- 0

#563-378 = 185 (fish families in fishTree but not in fish base)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","present_species_count_percentage")

#make a new dataframe by extracting family and missing percentage columns
df_merged_4 <- df_merged_2[,c("family","present_species_count_percentage")]

#Adding 185 families with 100% missing species percentage to 362 families, 
#the percentages are zeros here for present species
new_present <- rbind(df_merged_4,setdif_FishBase)
class(new_present)
#check the bottom of the dataset
tail(new_present)

attach(new)
#Histogram
hist(missing_percentage,main = "Histogram of percent missing species per family", xlab = "Missing species per family (%)")

#Statistics
sd(missing_percentage)
mean(missing_percentage)
median(missing_percentage)
var(missing_percentage)
length(missing_percentage)

detach(new)
attach(new_present)
#Histogram
hist(present_species_count_percentage ,main = "Histogram of percent present species per family", xlab = "Present species per family (%)")
abline(v = median(present_species_count_percentage),                     # Add line for median
       col = "red",
       lty=2,
       lwd = 2)
text(x = median(present_species_count_percentage) * 5,                 # Add text for median
     y = median(present_species_count_percentage) * 15,
     paste("Median =", median(present_species_count_percentage)),
     col = "red",
     cex = 1.2)
#Statistics
sd(present_species_count_percentage)
mean(present_species_count_percentage)
median(present_species_count_percentage)
var(present_species_count_percentage)
length(present_species_count_percentage)

#Checking the best distribution which is closest to the real dataset
set.seed(123)
x <- rgamma(n=547, shape= 0.2,  rate = 1/25.50842)
gamma <- fitdistr(x, "gamma")
summary(gamma)
hist(x,main= "Gamma Distribution")
gamma$loglik
AIC(gamma)

#To get the best shape value check the likelihood value
set.seed(123)
x <- rgamma(n=547, shape= 0.2,  rate = 1/25.50842)
gamma <- fitdistr(x, "gamma")
gamma$loglik


set.seed(123)
x4 <- rnegbin(547, mu = 25.50842, theta = 1)
neg_bino <- fitdistr(x4, "Negative Binomial")
neg_bino$loglik
AIC(neg_bino)
summary(neg_bino)
hist(x4, main= "Negative binomial Distribution")

set.seed(123)
x3 <- rweibull(547, shape = 1, scale = 25.50842)
wbul <- fitdistr(x3, "weibull")
wbul$loglik
AIC(wbul)
hist(x3,main= "Weibull Distribution")
summary(wbul,main= "Distribution of Weibull")

x5 <- rexp(n=547, rate=1/25.50842)
exp <- fitdistr(x5,"exponential")
hist(x5,main= "Exponential distribution")
exp$loglik
AIC(exp)
summary(exp)

#AIC comparison of models
AIC(neg_bino,wbul,exp,gamma)

#I'm using gamma to simulate missingness. As it's the best model for that according to the AIC comparison
#random gamma distribution (rate=(1/mean))

set.seed(1002)
#get 10 sets of 547 values using gamma distribution
sample_10_gamma <- data.frame(replicate(10,rgamma(n=547, shape =0.2, rate = 1/25.50842)))
head(sample_10_gamma)

#get 10 sets of 547 values using exponential distribution
sample_10_exponential <- data.frame(replicate(10,rexp(n=547, rate = 1/25.50842)))
head(sample_10_exponential)

#change values over 100 to 100
sample_10_gamma[sample_10_gamma > 100] <- 100
class(sample_10_gamma)

#change values over 100 to 100
sample_10_exponential[sample_10_exponential > 100] <- 100
class(sample_10_exponential)

# #In the stratified sampling we had 378 families. So here I'm gonna use 378 random values based on the model
# sample_10_gamma_1 <- data.frame(sample(nrow(sample_10_gamma),378))
# sample_10_gamma <- as.data.frame(sample_10_gamma)

#now using the 10 sets of random numbers do the distribution
for (col in 2:ncol(sample_10_gamma)) {
  hist((sample_10_gamma[,col]),main = "Gamma distribution for present species", xlab = "Gamma random numbers")
}
class(sample_10_gamma)

#now using the 10 sets of random numbers do the distribution
for (col in 2:ncol(sample_10_exponential)) {
  hist((sample_10_exponential[,col]),main = "Exponential distribution for present species", xlab = "Exponential random numbers")
}
class(sample_10_exponential)

#Bind the dataframe with random numbers (sample_10_present) to the dataframe consists of family name and species count
binded_df <- cbind(new_present,sample_10_gamma)
tail(binded_df)

#Bind the dataframe with random numbers (sample_10_present) to the dataframe consists of family name and species count
binded_df_expo <- cbind(new_present,sample_10_exponential)
tail(binded_df_expo)

#get species count to the dataframe (Here,the families with zero species are left out)
merged_df <- merge(df_tt,binded_df,by="family")

class(merged_df)
detach(new_present)
attach(merged_df)
head(merged_df)
#Calculate the sample species to extract from the dataset
#merged_df1<- merged_df %>% mutate(sampled_species = (round(species_count*X1/100)))

#merged_df1<- merged_df %>% mutate(sampled_species = (ceiling(species_count*X1/100)))
#head(merged_df1)
merged_df1_present <- round(merged_df[,2]*merged_df[,4:13]/100)
head(merged_df1_present)
merged_df1_missing <- merged_df[,2]- merged_df1_present[,1:10]
head(merged_df1_missing)
tail(merged_df1_missing)

miss_wt_families <- cbind(merged_df$family,merged_df1_missing)
head(miss_wt_families)
#first we need present species to build trees
present_wt_families <- cbind(merged_df$family,merged_df1_present)
head(present_wt_families)

#Rename column names
colnames(present_wt_families) <- c("family","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")

#checking
head(present_wt_families)
class(present_wt_families)
# #get the first biased sample with family name
# present_wt_families_1 <- present_wt_families[,1:2]

#create a function to split the large data frame into small data frames

split <- function(column){
  df_split <- present_wt_families[,c(1,column)]
  return(df_split)
}

library(foreach)
# #A list of columns
# col_vector <- c(2:11)
# 
# #we can use foreach loop to iterate through each data frame (There were 40 data frames)
# split_col <- foreach(i=1:length(col_vector)) %do% present_wt_families[,c(1,col_vector[[i]])]
#   
# head(split_col)


#call the split function to generate 10 data frames splitting a large data frame
df_01 <- split(2)
df_02 <- split(3)
df_03 <- split(4)
df_04 <- split(5)
df_05 <- split(6)
df_06 <- split(7)
df_07 <- split(8)
df_08 <- split(9)
df_09 <- split(10)
df_10 <- split(11)

#checking
head(df_10)
class(df_10)
#function to remove families with zero species
fam_no_zero <- function(df_new,C){
  new_df<- df_new[df_new$C >0,]
  return(new_df)
}


# call the fam_no_zero funtion to generate backbone trees that represent present species according to biased realistic samples
d01 <- fam_no_zero(df_01,"C1")
d02 <- fam_no_zero(df_02,"C2")
d03 <- fam_no_zero(df_03,"C3")
d04 <- fam_no_zero(df_04,"C4")
d05 <- fam_no_zero(df_05,"C5")
d06 <- fam_no_zero(df_06,"C6")
d07 <- fam_no_zero(df_07,"C7")
d08 <- fam_no_zero(df_08,"C8")
d09 <- fam_no_zero(df_09,"C9")
d10 <- fam_no_zero(df_10,"C10")
class(d10)

#pull function of the dplyr package to convert a column of a data frame into a vector
#C1 <- pull(present_wt_families_1,C1)
C1 <- pull(d01,C1)
C2 <- pull(d02,C2)
C3 <- pull(d03,C3)
C4 <- pull(d04,C4)
C5 <- pull(d05,C5)
C6 <- pull(d06,C6)
C7 <- pull(d07,C7)
C8 <- pull(d08,C8)
C9 <- pull(d09,C9)
C10 <- pull(d10,C10)
class(C1)
# sample_sizes <- data.frame(
#   families = unique(present_wt_families_1$family),
#   n_to_sample = C1
# )
# 
# sample_size <- function (df,C){
#   sample_df <- data.frame(families = unique(df$family), n_to_sample = C)
#  return(sample_df) 
# }
# 
# df_1s <- sample_size(d01,C1)

#join two data frames together. So now sequences also in the same data frame
join_d01 <- inner_join(d01,df_fish,by="family") 
class(join_d01)
tail(join_d01)
join_d02 <- inner_join(d02,df_fish,by="family") 
join_d03 <- inner_join(d03,df_fish,by="family") 
join_d04 <- inner_join(d04,df_fish,by="family") 
join_d05 <- inner_join(d05,df_fish,by="family") 
join_d06 <- inner_join(d06,df_fish,by="family") 
join_d07 <- inner_join(d07,df_fish,by="family") 
join_d08 <- inner_join(d08,df_fish,by="family") 
join_d09 <- inner_join(d09,df_fish,by="family") 
join_d10 <- inner_join(d10,df_fish,by="family") 
# join <- function(df){
#   df_seq <- inner_join(df,df_fish,by="family")
#   return(df_seq)
# }

# group <- function(df_j,C){
#   df_join <- df_j %>% 
#     group_by(family) %>% 
#     sample_n(C)
#   return(df_join)
#   
# }
# 
# j1 <- join (d01)
# g <- group(j1,C1)

#get different samples based on family. This contains the sequences as well
df_g1 <- join_d01 %>% 
  group_by(family) %>% 
  sample_n(C1)

class(df_g1)

df_g2 <- join_d02 %>% 
  group_by(family) %>% 
  sample_n(C2)

df_g3 <- join_d03 %>% 
  group_by(family) %>% 
  sample_n(C3)

df_g4 <- join_d04 %>% 
  group_by(family) %>% 
  sample_n(C4)
df_g5 <- join_d05 %>% 
  group_by(family) %>% 
  sample_n(C5)

df_g6 <- join_d06 %>% 
  group_by(family) %>% 
  sample_n(C6)
df_g7 <- join_d07 %>% 
  group_by(family) %>% 
  sample_n(C7)

df_g8 <- join_d08 %>% 
  group_by(family) %>% 
  sample_n(C8)
df_g9 <- join_d09 %>% 
  group_by(family) %>% 
  sample_n(C9)

df_g10 <- join_d10 %>% 
  group_by(family) %>% 
  sample_n(C10)

#The data type should be dataf.frame to use dat2phylip (). otherwise the result give only one sequence in phylip output file. Soconvert to dataframe
df_g1 <- as.data.frame(df_g1)
class(df_g1)
df_g2 <- as.data.frame(df_g2)
class(df_g2)
df_g3 <- as.data.frame(df_g3)
df_g4 <- as.data.frame(df_g4)
df_g5 <- as.data.frame(df_g5)
df_g6 <- as.data.frame(df_g6)
df_g7 <- as.data.frame(df_g7)
df_g8 <- as.data.frame(df_g8)
df_g9 <- as.data.frame(df_g9)
df_g10 <- as.data.frame(df_g10)
#create a function to subset data frames 

df_sub<- function(df){
  df_subset <- df[,c(3:4)]
  return(df_subset)
}

df1 <- df_sub(df_g1)
class(df1)
df2 <- df_sub(df_g2)
class(df2)
df3 <- df_sub(df_g3)
df4 <- df_sub(df_g4)
df5 <- df_sub(df_g5)
df6 <- df_sub(df_g6)
df7 <- df_sub(df_g7)
df8 <- df_sub(df_g8)
df9 <- df_sub(df_g9)
df10 <- df_sub(df_g10)

library(phylotools)
#Use the function to generate phylip files from the dataframe
dff <- function(xx){
  return(dat2phylip(xx,outfile = "bias.phy"))
}

dff(df10)

tail(group)
df_fish_present <- inner_join(df_fish,present_wt_families)
tail(df_fish_present)

df_fish_present <- merge(df_fish,present_wt_families,by="family",all=T)
head(df_fish_present)
tail(df_fish_present)

wot <- rep(present_wt_families$family,present_wt_families$C1)
tail(wot)



