#### PART 0: INFORMATION --------
# Written by Alice Wang, ID 1080579, in fufilment of BINF*6210 Assignment 2 on the use of ND1 sequence data on frog (Order: Anura)

#### PART 1: PACKAGES ------
#Installing (if not previously installed, remove #) and loading all relevant libraries.
#install.packages("tidyverse")
library(tidyverse)
#install.packages("randomForest")
library(randomForest)
#install.packages("rentrez")
library(rentrez)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("seqinr")
library(seqinr)
#install.packages("caTools")
library(caTools)
#install.packages("caret")
library(caret)
#install.packages("MLmetrics")
library(MLmetrics)

##FUNCTIONS 
#pie of genus' to be represented in classifier after various steps of filtering 
pie_of_genus <- function(dataframename) {
  pie(dataframename$count, labels = dataframename$Genus_Name, main = "Relative proportions of genus'", col = rainbow(length(dataframename$Genus_Name)))
}
#making a plot of the importance of each variable in the predictor model 
make_importance_plot <- function(RF_model) {
  plot(RF_model$importance,
       main = "Importance",
       xlab = "k-mer Predictor",
       ylab = "Relative Importance (Mean Decimal in Gini Coefficient)",
       pch = 1,
       cex = 0.6,
       xaxt = "n")
  axis(1, at = (1:(length(RF_model$importance))), labels = rownames(RF_model$importance),tck = 0, cex.axis = 0.6, las = 2)
}

#### PART 2: RETRIVING INFORMATION ----- 
#ND1 is positioned within the Anura mitochondrial genome from position 2935 - 3897 which corresponds to a gene length of 963bp as per DOI: 10.7717/peerj.7532

#Retrive ND1 sequences of Anura that are at least 60% of the gene length (estimated 600bp). Allow sequences that are slightly longer than the gene (up to 1100bp).Set retmax to 250 to increase from the default 20. Initial search to find max n hits.
ND1_Ann <- entrez_search(db = "nuccore", term = "Anura[ORGN] AND ND1[GENE] AND 600:1100[SLEN]", retmax = 250)

#Make working dataset called ND1_Ann, and set retmax to max number of hits :)
ND1_Ann <- entrez_search(db = "nuccore", term = "Anura[ORGN] AND ND1[GENE] AND 600:1200[SLEN]", retmax = ND1_Ann$count, use_history = TRUE)

#Query for summaries for the hits in batches not to overwhelm the average laptop bc hits as of 10.25.23 is >4000 
sz_batch <- 100 # take 100 entries at a time. sz = size of batch

# split entries into batches of up to 100 
ids_batch <- split(ND1_Ann$ids, ceiling(seq_along(ND1_Ann$ids)/sz_batch))

#establish a vector to contain IDs after batching is complete !
Ann_ids <- vector(mode = "character", 0)

#loop through each batch to identify and store the IDs and append it to the Ann_ids vector
j = 1
while (j <= length(ids_batch)) {
  ND1_Ann_working <- entrez_summary(db = "nuccore", id = ids_batch[[j]])
  ids_working <- extract_from_esummary(ND1_Ann_working, "uid")
  Ann_ids <- c(Ann_ids, ids_working)
  j = j + 1
}

# separate entries into groups of100
ids_batch <- split(Ann_ids, ceiling(seq_along(Ann_ids)/sz_batch))

#create vector to contain sequences following batchin'
Ann_seq <- vector(mode = "character", 0) 

#loop through each chunk of IDs to retrieve sequence data for each and store it into the Ann_Seq vector
j = 1
while (j<= length(ids_batch)) {
  seq_working <- entrez_fetch(db = "nuccore", id = ids_batch[[j]], rettype = "fasta")
  Ann_seq <- c(Ann_seq, seq_working)
  j = j+1
}

#make file w sequences and make the delimitter a new line instead of default
write(Ann_seq, "ann_seq.fasta", sep = "\n")


#### PART THREE: QUALITY CHECKS AND FILTERING --------
#create a dataframe of sequences in string sets (sg set)
sg_set <- readDNAStringSet("ann_seq.fasta")
df_ND1_Ann <- data.frame(ND1_title = names(sg_set), ND1_Sequence = paste(sg_set))

#count length of each sequence by looping through all seq and adding them to numeric vector containing all seq lens 
seq_len <- vector(mode="numeric",0)

j = 1
while (j<= length(df_ND1_Ann$ND1_Sequence)) {
  seq_len <- c(seq_len, nchar(df_ND1_Ann$ND1_Sequence[j]))
  j = j+1
}

#visualize distribution of seq lengths 
hist (seq_len,
      main = "Distribution of Sequence Lengths",
      breaks = 10,
      xlab = "Sequence length (nt)",
      ylab = "Number of sequences",
      col = "pink",
      ylim = c(0,1200),
      las = 2)

#Obtain the genus under the order Ann which each sequence belongs to
#create new column to store genus name 
df_ND1_Ann$Genus_Name <- word(df_ND1_Ann$ND1_title, 2L)
#rearrange this dataframe
df_ND1_Ann <- df_ND1_Ann[, c("ND1_title", "Genus_Name", "ND1_Sequence")]

#inquiry into how many occurrences of each genus there are
n_genuss <- df_ND1_Ann %>%
  group_by(Genus_Name) %>%
  summarise(count=n())
View(n_genuss)
#table copy to move to doc/excel: write.table(n_genuss, "clipboard", sep="\t", row.names=FALSE)

#remove genera which are not represented well. Keep only genera that occur over 100 times.
gen_tokeep <- subset (n_genuss, count > 100)
df_filtered_genus <- df_ND1_Ann[df_ND1_Ann$Genus_Name %in% gen_tokeep$Genus_Name, ]
pie_of_genus(gen_tokeep)

#note: class imbalance ! (uhoh... address this later). Several not even shown on the chart because so small...

#remove sequence gaps in cases which the sequencing wasn't of high quality and there are a lot of undetermined nucleotides (N = unknown nucleotide). Remove N > 1% of the sequence
df_filtered_gap <- df_filtered_genus %>%
  mutate(ND1_Seq_2 = str_remove(ND1_Sequence, "^[-N]+")) %>%
  mutate(ND1_Seq_2 = str_remove(ND1_Seq_2, "[-N]+$")) %>%
  mutate(ND1_Seq_2= str_remove_all(ND1_Seq_2, "-+")) %>%
  filter(str_count(ND1_Seq_2, "N") <= (0.01 * str_count(ND1_Sequence)))

#to compare sequences of similar length, use only two quartiles of sequence length. 
#script used for mid quartl
MID_lowlim <- quantile(nchar(df_filtered_gap$ND1_Seq_2), probs = 0.25, na.rm = TRUE)
MID_uprlim <- quantile(nchar(df_filtered_gap$ND1_Seq_2), probs = 0.75, na.rm = TRUE)
MID_df_filtered_len <- df_filtered_gap %>%
  filter((str_count(ND1_Seq_2) >= MID_lowlim & str_count(ND1_Seq_2) <= MID_uprlim)) 
#create dataframe with filtered sequence data and genus identifications 
MID_df_gen_seq <- as.data.frame(MID_df_filtered_len [,c("Genus_Name", "ND1_Seq_2")])

MID_n2_genuss <- MID_df_gen_seq %>%
  group_by(Genus_Name) %>%
  summarise(count=n())
#remove underrepresented classes
MID_n2_genuss <- subset(MID_n2_genuss, count > 50)
pie_of_genus(MID_n2_genuss)
##appears that the class imbalance is even worse...

###NOW TRY WITH UPPER TWO Quartiles
lowlim <- quantile(nchar(df_filtered_gap$ND1_Seq_2), probs = 0.5, na.rm = TRUE)
df_filtered_len <- df_filtered_gap %>%
  filter((str_count(ND1_Seq_2) >= lowlim)) 
#create dataframe with filtered sequence data and genus identifications 
df_gen_seq <- as.data.frame(df_filtered_len [,c("Genus_Name", "ND1_Seq_2")])

n2_genuss <- df_gen_seq %>%
  group_by(Genus_Name) %>%
  summarise(count=n())
#remove severely underrepresented classes
n2_genuss <- subset(n2_genuss, count > 50)
df_gen_seq <- df_gen_seq[df_gen_seq$Genus_Name %in% n2_genuss$Genus_Name, ]

pie_of_genus(n2_genuss)
view(n2_genuss)
#ah, much better :) but we will still have to count for class imbalance later!
rm(Ann_ids, Ann_seq, df_filtered_gap, df_filtered_genus, df_filtered_len, df_ND1_Ann, gen_tokeep, ids_batch, ids_working, lowlim, MID_df_filtered_len, MID_lowlim, MID_n2_genuss, MID_uprlim, n_genuss, n2_genuss, ND1_Ann, ND1_Ann_working, seq_len, seq_working, sg_set, sz_batch)

#### PART FOUR: ADDING FEATURES --------------
#take the sequences and convert into  DNAStringSet for Biostring use
df_gen_seq$ND1_Seq_2 <- DNAStringSet(df_gen_seq$ND1_Seq_2)

#Adding each nucleotide's frequency 
df_gen_seq_plnt <- cbind(df_gen_seq, as.data.frame(letterFrequency(df_gen_seq$ND1_Seq_2, letters = c("A", "T", "C", "G"))))

#adding dinuctide and trinucleotide frequencies 
df_gen_seq_aftdi <- cbind(df_gen_seq_plnt, as.data.frame(dinucleotideFrequency(df_gen_seq_plnt$ND1_Seq_2, as.prob = TRUE)))

df_gen_seq_aftri <-cbind(df_gen_seq_aftdi, as.data.frame(trinucleotideFrequency(df_gen_seq_aftdi$ND1_Seq_2, as.prob = TRUE)))
rm(df_gen_seq_aftdi)
df <- df_gen_seq_aftri[-c(2)]
df <- data.frame(df)
#to check if the dataframe looks correct and ncol should be 85
#view(df)
#ncol(df)

#to address the class imbalance issue... originally considered upSample from the carat package. however, this is best when there are only 2 classes. Left in script to show consideration, but to explain why it was not used.
#up_trn <- upSample(x = trn[-1], y=trn$Genus_Name)
#SMOTE was also considered but not validated in biological contexts. 

#now to address the minority class...  we will up sample (randomly) to get to closer to the size of the other classses. Will also randomly (down) sample the overrepresentated class.

#can be altered - create sample size to be added
sample_sz <- 80
#dataset to draw our upsampling from - from the tiny class 
prtdata <- subset(df, df$Genus_Name == "Pristimantis")
#creating indexing variable to make random sampling numbers to choose which rows to upsample!
tot_rws <- nrow(prtdata)
spd_inx <- sample (1:tot_rws, size = sample_sz, replace = TRUE)
up_prt <- prtdata[spd_inx, ]
df <- rbind(df, up_prt)

#for the down-sample - we will max out at 250 to match the other classes
sample_sz <- 250
#dataset to draw our sampling from
scidata <- subset(df, df$Genus_Name == "Scinax")
#creating indexing for rows and randomly sampling rows to upsample!
tot_rws <- nrow(scidata)
spd_inx <- sample (1:tot_rws, size = sample_sz, replace = TRUE)
down_prt <- scidata[spd_inx, ]


#removing all Scinax from df
df <- subset(df, df$Genus_Name != "Scinax")
#adding back the sampled scinax
df <- rbind(df, down_prt)

#final genus composition to use 
end_genuss <- df %>%
  group_by(Genus_Name) %>%
  summarise(count=n())
pie_of_genus(end_genuss)
#view(end_genuss)

#Splitting our dataframe into two subsets - 80% of the data will be used to train our model. Then 15% will be reintroduced to the model (without the results) to test how well our model works
#Setting seed to my student ID so this can be replicated.
set.seed(1080579)
df$Genus_Name <- as.factor(df$Genus_Name)
train_index <- sample(1:nrow(df),0.8*nrow(df))
trn <- df[train_index,]
tst <- df[-train_index,]

#the 3-mer columns to be used as predictors
predictor_idx_tri <- c(22:85)
#selecting the 80% predictors to be used for training
predictor_tri <- df[,predictor_idx_tri]

#column 1 is genus name. #column 2-6 is single nt freq. #col 6-21 is is dint. #col 22-85 is trinuct.
#k-mer rows calculated by:
#START: sum of 4^(n-1) + 1 for n = (k to 0) where k is the number of nucleotides in the k-mer.
#END: sum of 4^(n) for n = (k to 0) where k is the number of nucleotides in the k-mer.
#ex: for 3-mers (trinucleotides), rows are... #START: (4^2 + 4^1 + 4^0 + 1) = 22. END: (4^3 + 4^2 + 4^1 + 4^0) = 85. 

#creating a random forest classifier for 3-mers and creating 25 trees.
RF_clsfr <- randomForest(x = predictor_tri,
                         y = df$Genus_Name,
                         ntree =25)
RF_clsfr
#testing the 3-mer classifier
y_pd <- predict(RF_clsfr, newdata = tst[-1])
mx_confusion <- confusionMatrix(y_pd, reference = tst[,1])
mx_confusion
#wow, perfect!

#plotting the relative importance of each of the 3-mers
#write.table(RF_clsfr$importance, "clipboard", sep="\t", row.names=TRUE)
view(RF_clsfr$importance)
make_importance_plot(RF_clsfr)

#PART 6: SECOND MODEL, only 2mer ----------------------------------
#using the correct rows in the dataframe that use 2-mers (dinucleotides)
predictor_idx_di <- c(6:21)
predictor_di <- df[,predictor_idx_di]
#set classifier
RF_clsfr2 <- randomForest(x = predictor_di,
                          y = df$Genus_Name,
                          ntree =25)
RF_clsfr2

#test!
y_pd2 <- predict(RF_clsfr2, newdata = tst[-1])
mx_confusion2 <- confusionMatrix(data=y_pd2, reference = tst[,1])
mx_confusion2
#wow, perfect again!

#view importance of each dinucleotide
#view(RF_clsfr2$importance)
make_importance_plot(RF_clsfr2)
#write.table(RF_clsfr2$importance, "clipboard", sep="\t", row.names=TRUE)


#### PART 7: SINGLE NT CLASSIFIER ?!
#select only single nucleotide frequencies 
predictor_idx_nt <- c(2:5)
predictor_nt <- df[,predictor_idx_nt]
#build classifier
RF_clsfr3 <- randomForest(x = predictor_nt,
                          y = df$Genus_Name,
                          ntree =25)
#test it 
y_pd3 <- predict(RF_clsfr3, newdata = tst[-1])
mx_confusion3 <- confusionMatrix(data=y_pd3, reference = tst[,1])
mx_confusion3
#perfect :D

#view importance
view(RF_clsfr3$importance)
#write.table(RF_clsfr3$importance, "clipboard", sep="\t", row.names=TRUE)
make_importance_plot(RF_clsfr3)

#############PART 8: COMPARING THEM ALL-----------

#compiling the metrics of each classifier in order of 3-mer, 2-mer, nt frequency (the order in which they were built)
metrics_lst <- list(
  RF_clsfr = mx_confusion$overall,
  RF_clsfr2 = mx_confusion2$overall,
  RF_clsfr3 = mx_confusion3$overall
)

#create a dataframe to show each of these models' metriccs side by side for easy comparision 
df_metrics <- do.call(cbind, lapply(names(metrics_lst), function(model) {
  data.frame(Model = model, metrics_lst[[model]])
}))

#editing names
df_metrics_sum <- df_metrics[,c(2,4,6)]
colnames(df_metrics_sum)<- c("3-mer", "2-mer", "nt")

#let's take a look!
view(df_metrics_sum)
#write.table(df_metrics_sum, "clipboard", sep="\t", row.names=TRUE)

#looks like they all perform amazingly. 
