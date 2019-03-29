#############################################################
# Ref to the ARTICLE
# 
# Pilar Morera Margarit & Davide Bulgarelli
# Pilar.MoreraMargarit@hutton.ac.uk
# d.bulgarelli@dundee.ac.uk
# 
# script to reproduce calculations and figures presented in the manuscript

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("PMCMR")
#install.packages("wesanderson")
#install.packages("dendextend") #you need to say no when it asks you

#biocLite ("exactRankTests")
#install package ANCOM from zip file

#required packages 
library("phyloseq")
library("PMCMR")
library("vegan")
library ("ape")
library("dendextend")
library("ggplot2")

#Install ancom function for Mac
source ("ancom_functions.R")
source("plot_ancom.R")
#Install ancom function for Windows
#library("ancom.R")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory
#Pilar cpu
setwd("/Users/Pilar/R sequencing analysis")

#############################################################
#Import the count matrix and the desing file
#############################################################
#OTU table generated using QIIME 1.9.0. 
dat_info <- read.delim("JH09_PM_otu_table_nc2_checked.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#Inspect the file 
dim(dat_info)
colnames(dat_info)

#Extract the total number of reads clustered at OTU 97% identiy 
#(the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:44]))
OTU_97_reads

#Total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:44]))
OTU_97_reads_sum

#Design file
design <- read.delim("Map_JH09_PM.txt", sep = "\t", header=TRUE, row.names=1)
design

#Remove chloroplast and mitochondria OTUs from the original dataset
#Inspect the first rows of the table to get insights into the column ConsensusLineage
dat_info[1:5, 45]

#Subset the original OTU table and create a new one with:
#OTUs assigned to chloroplast
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$ConsensusLineage), ]
dim(Chloroplast)

#and OTUs assigned to mitochondria
mitochondria <- dat_info[grepl("mitochondria", dat_info$ConsensusLineage), ]
dim(mitochondria)

#Set a difference between the row names of the the three datasets: 
#This information will be used to filter out Plant and Host derived OTUs from the OTU table
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))

#Inspect for potential contaminants 
#(data based on the results from Laura Pietrangello)
contaminants <- read.delim("contaminant_OTUs_30116224.txt", header = T, row.names = 1)
rownames(contaminants)

#Create a dataset without the conatminants 
#We use noPlants so the result will also be without chloroplasts and mitochondria sequences
#name noPC means no Plants and Contaminants
noPC <- setdiff(noPlants, rownames(contaminants))

#Inspect the results
length(rownames(dat_info))
length(noPlants)
length(noPC)

#Save the OTUids list generated at line 100 
#This will be used to usbset the OTU table and generate taxa-tables in QIIME
#write(noPC, "JH09_PMM_noPC_OTUs_id.txt")

#Generate a new OTU table which will be devoid of Chloroplast and Mitochondria OTUs
#After revision of the manuscript we decided to remove contaminant OTUs from the dataset
#So, to not change all the names in the script, 
#we will keep the name noPlants instead of changing to noPC but 
#for now on noPlant will be the dataset deprived of chlroplast, mitochondria and contaminant sequences
dat_info_noPlants <- dat_info[noPC, ]

#Create a new count matrix without OTUs assigned to Choloplast, Mitochondria and Contaminants
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPC, rownames(design)]
dim(dat_count_noplants)
dat_count_noplants[1:50, ]
#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 45])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#Save the above file and in excel we will create a new taxa table where each column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH09_PMM_dat_tax_noPC.txt", sep="\t")

#Check the effect of mitochondria, chloroplast and contaminants depletion on the new OTU table
#Dataset with chloroplast, mitochondria and contaminant sequences
dim(dat_count)

#Dataset without chloroplast, mitochondria and contaminant sequences
dim(dat_count_noplants)

#Total number of reads without chloroplast, mitochondria and contaminant sequences per sample
OTU_97_reads_noPlants <- colSums(dat_count_noplants)
OTU_97_reads_noPlants

#Now sort the samples in increasing order of reads
sort(OTU_97_reads_noPlants)

#Now calculate the total number of reads
OTU_97_reads_noPlants_sum <- sum(OTU_97_reads_noPlants)
OTU_97_reads_noPlants_sum 

#Define the proportion of non-chloropolast, mitochondria and contaminant reads in the original dataset 
useful_reads <- (OTU_97_reads_noPlants_sum/OTU_97_reads_sum)*100
useful_reads #91% of the reads were kept after filtering

#Create a dataset to visualise the proportion of reads per sample 
#before and after removing OTUs assigned to chloroplast, mitochondria and contaminants
OTU_97_reads_noPlants <- as.data.frame(OTU_97_reads_noPlants)
OTU_97_microbial_reads_proportion <- as.data.frame(colSums(dat_count_noplants)/colSums(dat_count))*100

#Rename the columns in the generated datasets
colnames(OTU_97_reads_noPlants) <- c("reads")
colnames(OTU_97_microbial_reads_proportion) <- c("Microbial_OTUs_reads")

#Combine these datasets with the design file
design_info_2 <- cbind(design, OTU_97_reads_noPlants, OTU_97_microbial_reads_proportion)

#Calculate the max, min and mean number of reads for the dataset
mean(design_info_2$reads)
max(design_info_2$reads)
min(design_info_2$reads)

######################################################################################
                          #Reads distribution
####################################################################################
#Data visualisation: Number of reads per sample
#reoreder the levels of a given factor for graphical representation
design_info_2$Location <- ordered(design_info_2$Location, levels=c("Stafford1", "Stafford2", "Shifnal", "Woore", "Invergowrie1", "Invergowrie2")) 

#Boxplot: Number of reads per population
with(design_info_2, boxplot(Microbial_OTUs_reads ~ Location, xlab = "Populations", ylab = "Sequencing Reads Proportion",   main = "Reads assigned to Microbial OTUs"))

#Test the normality of the observed data
shapiro.test(design_info_2$Microbial_OTUs_reads) #Data is not normally distributed

#Use a non parametric anova 
kruskal.test(Microbial_OTUs_reads ~ Location, data = design_info_2)
posthoc.kruskal.dunn.test (x=design_info_2$Microbial_OTUs_reads, g=design_info_2$Location, p.adjust.method="BH")
#############################################################
            #Visualize the taxa distribution  
############################################################

#***************************** Genus level ******************************
#Import average value genus level
dat_info_taxa_Genus <- read.delim("JH09_PM_otu_table_nc2_checked_manually_annotated_noPC_3_L6.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Genus)
dim(dat_info_taxa_Genus)

#Transform the data in % for visualisation
dat_norm_Genus <- dat_info_taxa_Genus * 100
dat_norm_Genus[1:5, ]
colSums(dat_norm_Genus)

#Determine the average % of reads for each genus
Genus_mean_sorted <- dat_norm_Genus[(order(-rowSums(dat_norm_Genus))), ]
Genus_mean_sorted

#Identify the genus whose average abundance is above 1% and their aggregated relative abundance
Genus_mean_topRank <- Genus_mean_sorted[rownames(Genus_mean_sorted)[which(rowMeans(Genus_mean_sorted) > 1)], ]
dim(Genus_mean_topRank)
colSums(Genus_mean_topRank)
rowMeans(Genus_mean_topRank)

#Transform in matrix for plotting purposes
Genus_mean_topRank <- as.matrix(Genus_mean_topRank)
colnames(Genus_mean_topRank)

#re-arrange the order for plotting purposes
col.order <- c("Strawberrry1I", "Strawberrry1II", "Strawberrry1III", "Strawberrry1IV", "Strawberrry1V", "Strawberrry1VI", "Strawberrry1VII", "Strawberrry1VIII",
               "Strawberrry2I", "Strawberrry2II", "Strawberrry2III", "Strawberrry2IV", 
               "TP07I","TP07II", "TP07III", "TP07IV", "TP07V", "TP07VI", "TP07VII", "TP07VIII",
               "TP08I","TP08II", "TP08III", "TP08IV", "TP08V", "TP08VI", "TP08VII", "TP08VIII",
               "PMM01I", "PMM01II", "PMM01III", "PMM01IV", "PMM01V", "PMM01VI", "PMM01VII", "PMM01VIII",
               "Strawberrry3I", "Strawberrry3II", "Strawberrry3III", "Strawberrry3IV", "Strawberrry3V", "Strawberrry3VI", "Strawberrry3VII", "Strawberrry3VIII")
Genus_mean_topRank <- Genus_mean_topRank[,col.order]

# Stacked Bar Plot with Colors and Legend
barplot(Genus_mean_topRank, 
        main="Genus Distribution",
        xlab="Populations", 
        ylab = "Reads %", 
        ylim = c(0,100), 
        col=c("coral", "darkorchid1", "blue"), 
        beside=FALSE,   
        legend = rownames(Genus_mean_topRank),
        args.legend = list(x = "topright", bty = "n", inset=c(-0.3, 0)))

#Due to size limits the legend covers part of the graph, omit the legend and add it a posteriori
jpeg(filename = "./figures/Barplot_genus_Nard.jpg", 
     res = 300,
     width = 224, 
     height = 174, 
     units = "mm")
par(family="Arial", cex.lab=1.2, 
    xpd= TRUE, cex=1.2, mgp=c(2.5,1,0), mar=c(4,4,1,0))
barplot(Genus_mean_topRank, 
        xlab="Populations", 
        ylab = "Reads %", 
        ylim = c(0,100), 
        col=c("black", "white", "grey"), 
        beside=FALSE,   
        bty = "n",
        axisnames = FALSE)
segments(0, -3, 9.5, -3, col="black")
text(5, -7, "St1", col = "black")
segments(9.8, -3, 14.4, -3, col="black")
text(12, -7, "St2", col = "black")
segments(14.7, -3, 24, -3, col="black")
text(19, -7, "Shf", col = "black")
segments(24.3, -3, 33.6, -3, col="black")
text(29, -7, "W", col = "black")
segments(33.9, -3, 43.2, -3, col="black")
text(38.6, -7, "I1", col = "black")
segments(43.5, -3, 53, -3, col="black")
text(48, -7, "I2", col = "black")

rect(xleft = 34, xright = 52, 
     ybottom = 0, ytop = 22, 
     lwd= 3, 
     col="white")

par(font="3")
text(36,4, "Candidatus Nardonella",
     adj = c(0,0), col = "black")
text(36,10, "Rickettsia",
     adj = c(0,0), col = "black")
par(font="1")
text(36,16, "Rickettsiaceae",
     adj = c(0,0), col = "black")
points(35,5, pch= 15, cex= 1.2)
points(35,11, pch= 15, cex= 1.2, col="grey")
points(35,17, pch= 0, cex= 1.2)
dev.off()
#************************** Candidatus Nardonella abundance *****************

#What is the exact name on the dat_count table?
#We checked opening the dat_count on excel 
#New.ReferenceOTU0 is the name

#Subset for Nardonella OTU
#Create a dataset with only Nardonella by subseting the count matrix from the row that is called New.ReferenceOTU0
dat_count_noplants_Nard <- as.data.frame(dat_count_noplants["New.ReferenceOTU0", ])
#Inspect that dataset
dat_count_noplants_Nard
#Transform into relative abundance dividing by the total number of reads and make the %
dat_count_noplants_Nard_RA <- (dat_count_noplants_Nard/colSums(dat_count_noplants))*100
colSums(dat_count_noplants)
#Inspect that dataset
dat_count_noplants_Nard_RA
#Inspect the names of the rows of the design (names of OTUs
#and the names of the columns of the count matrix. 
#They should be the same and in the same order
rownames(design)
colnames(dat_count_noplants_Nard)
#Just in case the order was not the same, 
#we impose the order of the rows of the design into the new count matrix for Nardonella
#These names will be the columns [rows, columns]
dat_count_noplants_Nard_RA <- dat_count_noplants_Nard_RA[ ,rownames(design)]
#We need to transpose the generated dataset 
#to have the insects as rows and the Nardonella counts as columns 
#for plotting purposes with the function 't'
dat_count_noplants_Nard_RA_T <- t(dat_count_noplants_Nard_RA)
#Inspect the transposed dataset
dat_count_noplants_Nard_RA_T
#Combine the design and the transposed count file
#We want to have the information of Location inside for plotting purposes
#with the function 'cbind'
dat_count_noplants_Nard_RA_T_info <- cbind(design, dat_count_noplants_Nard_RA_T) 
#Inspect that dataset
dat_count_noplants_Nard_RA_T_info
#Mean % of sequencing reads assigned to Nardonella
mean(dat_count_noplants_Nard_RA_T_info$New.ReferenceOTU0)

########################################################################################
                  #Genererate the phyloseq object
########################################################################################
#Data required: 
#dat_count_noplants; 
#design, 
#JH09_PMM_dat_tax_noPlants_ordered.txt, and 
#JH09_PMM.tre
#*************************************************************************************
#Generate a phyloseq object: this is a file that includes all the relevant information for the analysis: the OTU table, the taxonmy information, the phylogenetic tree and a mapping file 

#The phyloseq package and the phyloseq object is just gorgeous. Because it creates a file that recapitulates 
#a) The OTU Table counts. This is a n rows x m columns dataframe, where n represents the individual OTUs clustered at 97% sequence similarity and m the samples. The cells of the dataframe are "filled" with the sequencing reads (often reported as just reads) generated in MiSeq and assigned to that OTU in given a given sample
#b) The taxonomy information. For each of the n OTUs, a representative sequence (the most abudant in the cluster) has been used to query a database containing the taxonomy information for Prokaryotes and identify the bugs closely related to that given OTUs
#c) The mapping file (oftend reported as a design file or metadata): it says what the samples are, from where they come from and also it includes additional attributes that can be used in the analysis to correlate microbiota data with host traits (e.g., treatment, dry weight...)
#d) The phylogenetic tree. A phylogenetic tree is a graphical representation of the evolutionary relationships among organisms. In our partcular case we deal with a "gene tree". 
#A tree is composed by two elements: nodes, representing the feature of the organisms under study (in our case individual OTUs), and branches connecting nodes, defining the relationships among nodes (e.g., sequence diversity among the 16S rRNA gene sequences we are analysing using)

#Now the two big advantages of phyloseq are: 
#1) once created a phyloseq object, you can manipulate simoultaneously all the relevant components at once (e.g., removing samples, removing OTUs) something that otherwise will take you ages and quite possibly a lot of mistake
#2) the second advantage is that phyloseq is a wrapper that can be used directly for a variety of data visualisation and calculation (this will be illustrated below)

#a)The OTU Table counts
JH09_PMM_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)

#b)The taxonomy information
#Note that the file JH09_PMM_dat_tax_noPC_ordered.txt has been generated from the output of lines 347-348  
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
JH09_PMM_taxa_ordered <- read.delim ("JH09_PMM_dat_tax_noPC_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH09_PMM_taxa <- tax_table(as.matrix(JH09_PMM_taxa_ordered))
dim(JH09_PMM_taxa)

#c)The mapping file 
JH09_PMM_map <- sample_data(design)

#d)The phylogenetic tree: this is the phylogenetic tree generated by QIIME for open ref
JH09_PMM_tree <- read_tree("JH09_PMM.tre")

#Assign a root to the tree
#Check whether the tree is rooted
is.rooted(JH09_PMM_tree)

# Is there any unique class in the database?
unique_OTUs <- unique(JH09_PMM_taxa_ordered[,2])
unique_OTUs

#We could use a Bacteroidetes as an outgroup 
outgroup <- JH09_PMM_taxa_ordered[grepl("p__Bacteroidetes", JH09_PMM_taxa_ordered$Phylum), ]

#How many of these p__Bacteroidetes do we have? 
dim(outgroup)

#Identify the top abundant OTU of this group
sort(rowSums(dat_count_noplants[rownames(outgroup), ]))

#bit confusing but with 745 reads, OTU 4455618 is relatively abundant and part of Bacteroidetes, let's pick this OTU as an outgroup
newRoot = c("4455618")
JH09_PMM_tree <- root(JH09_PMM_tree,newRoot,resolve.root = TRUE)
#is it rooted now?
is.rooted(JH09_PMM_tree)# Yes it is!

#Merge the files and create the phyloseq object
JH09_PMM_data_phyloseq <- merge_phyloseq(JH09_PMM_OTU, JH09_PMM_taxa, JH09_PMM_map,  JH09_PMM_tree)

#Inspect the generated data
JH09_PMM_data_phyloseq #I still have the same amount of OTUs

#Number of reads in the phyloseq object
sum(colSums(otu_table(JH09_PMM_data_phyloseq)))

#Initial number of reads
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

#Abundance filtering: 
#Remove OTUs tallyig less than 5 reads in at least 10% of the samples (=1 Location)
JH09_PMM_data_phyloseq_2 = filter_taxa(JH09_PMM_data_phyloseq, function(x) sum(x > 5) > (0.1*length(x)), TRUE)
JH09_PMM_data_phyloseq_2

#Proportion of retained OTUs
(length(rownames(otu_table(JH09_PMM_data_phyloseq_2))))/(length(rownames(otu_table(JH09_PMM_data_phyloseq)))) * 100

#Proportion of retained reads
(sum(colSums(otu_table(JH09_PMM_data_phyloseq_2))))/(sum(colSums(otu_table(JH09_PMM_data_phyloseq)))) *100

#Number of reads after abundance filtering
sum(colSums(otu_table(JH09_PMM_data_phyloseq_2)))

########################################################################################
                 #Alphadiversity calculations
########################################################################################
#Data required: 
#design; 
#JH09_PMM_data_phyloseq_rare_table_counts2.txt

#Rarefy the dataset
#JH09_PMM_data_phyloseq_rare <- rarefy_even_depth(JH09_PMM_data_phyloseq_2, rngseed=TRUE)

#Extract and save the OTU table for reproducibility of the code
#JH09_PMM_data_phyloseq_rare_table <- as.data.frame(otu_table(JH09_PMM_data_phyloseq_rare))

#Inspect the generated file
#class(JH09_PMM_data_phyloseq_rare_table)
#dim(JH09_PMM_data_phyloseq_rare_table)

#Save the file for the reproducibility of the code
#write.table(JH09_PMM_data_phyloseq_rare_table, file="JH09_PMM_data_phyloseq_rare_table_counts_PC.txt", sep="\t")

#Import the rarefied OTU counts 
#(note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("JH09_PMM_data_phyloseq_rare_table_counts_PC2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#Inspect the generated file
dim(dat_count_rare)
colSums(dat_count_rare) #indeed it has been rarefied as all the samples have the smame number of reads

#Generate a new phyloseq object wich will contain only the rarefied counts and the design file 
#only these two pieces of information are required for alphadiversity calculation
JH09_PMM_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
JH09_PMM_map_rare <- sample_data(design)
JH09_PMM_data_rare_phyloseq <- merge_phyloseq(JH09_PMM_OTU_rare, JH09_PMM_map_rare)

#Inspect the generated file
JH09_PMM_data_rare_phyloseq
sample_sums(JH09_PMM_data_rare_phyloseq)

#Index calculations
JH09_PMM_alpha_rare <-  estimate_richness(JH09_PMM_data_rare_phyloseq, measures = c("Observed", "Shannon", "Chao1"))
JH09_PMM_alpha_rare

#Merge the datasets design file and alpha rare file
JH09_PMM_alpha_rare_info <- cbind(design, JH09_PMM_alpha_rare[rownames(design), ])

#Calculate the average
mean(JH09_PMM_alpha_rare_info$Observed)
mean(JH09_PMM_alpha_rare_info$Chao1)
mean(JH09_PMM_alpha_rare_info$Shannon)

#                     Test the normality of the observed data
#**********************************************************************************
#https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test
#Normality is assumed when the pvalue is higher than 0.05, i.e. null hypothesis is not rejected
shapiro.test(JH09_PMM_alpha_rare_info$Observed) #It is normal
shapiro.test(JH09_PMM_alpha_rare_info$Chao1) #It is normal
shapiro.test(JH09_PMM_alpha_rare_info$Shannon) #It is NOT normal

#ANOVA combined with Tukey test as a posthoc analyisis
#padj is the adjusted pvalue for multiple comparisons
#Check for multiple comparisons:
#https://stats.stackexchange.com/questions/253588/interpreting-tukeyhsd-output-in-r
#https://en.wikipedia.org/wiki/Multiple_comparisons_problem
#https://en.wikipedia.org/wiki/Post_hoc_analysis

#Observed OTUs
Observed_aov <- aov(Observed ~ Location, data = JH09_PMM_alpha_rare_info)
summary(Observed_aov)
TukeyHSD(Observed_aov) 

#Chao1 index
Chao1_aov <- aov(Chao1 ~ Location, data = JH09_PMM_alpha_rare_info)
summary(Chao1_aov)
TukeyHSD(Chao1_aov) 

#Kruskal Wallis test paired with Dunn test for posthoc analysis
#Shannon index
kruskal.test(Shannon ~ Location, data = JH09_PMM_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH09_PMM_alpha_rare_info$Shannon, g=JH09_PMM_alpha_rare_info$Location, p.adjust.method="BH")

#Reoreder the levels of a given factor for graphical representation
JH09_PMM_alpha_rare_info$Location <- ordered(JH09_PMM_alpha_rare_info$Location, levels=c("Stafford1", "Stafford2", "Shifnal", "Woore", "Invergowrie1", "Invergowrie2")) 

#Do the same graphs but put them together and with the stats labels
#boxplot observed OTUs
jpeg(filename = "./figures/Boxplot_a-diversity.jpg", 
     res = 300,
     width = 224, 
     height = 174, 
     units = "mm")

par(mfrow =c(3,1), 
    mar= c(0.5,4.4,0.5,.1), 
    oma = c(6,0,0,0) + 
      0.1, 
    bty= "l", 
    family= "Arial") 

with(JH09_PMM_alpha_rare_info, 
     boxplot(Observed ~ Location, 
             xaxt="n", 
             ylim = c(20, 75), 
             cex.axis = 1.7, 
             cex.lab= 1.7, 
             cex.names = 1.7, 
             ylab = "Observed OTUs"))
text(x= 0.37, y = 70, labels = "A", cex = 2)     
text(x= 1, y= 50, labels= "a,d", cex = 2)
text(x= 2, y= 46, labels= "a,d", cex = 2)
text(x= 3, y= 57, labels= "a,b", cex = 2)
text(x= 4, y= 63, labels= "b", cex = 2)
text(x= 5, y= 34, labels= "c", cex = 2) 
text(x= 6, y= 41, labels= "c,d", cex = 2)
#boxplot Chao1 OTUs
with(JH09_PMM_alpha_rare_info, 
     boxplot(Chao1 ~ Location, 
             xaxt="n", 
             ylim = c(20, 85), 
             cex.axis = 1.7, 
             cex.lab= 1.7, 
             cex.names = 1.7, 
             ylab = "Chao1 index"))
text(x= 0.37, y = 80, labels = "B", cex = 2)
text(x= 1, y= 60, labels= "a,d", cex = 2)
text(x= 2, y= 57, labels= "a,d", cex = 2)
text(x= 3, y= 66, labels= "a,b", cex = 2)
text(x= 4, y= 80, labels= "b", cex = 2)
text(x= 5, y= 45, labels= "c", cex = 2)
text(x= 6, y= 47, labels= "c,d", cex = 2)
#boxplot Shannon OTUs
with(JH09_PMM_alpha_rare_info, 
     boxplot(Shannon ~ Location, 
             ylim = c(0, 1.7), 
             cex.axis = 1.7, 
             cex.lab= 1.7, 
             cex.names = 1.7, 
             ylab = "Shannon index", 
             xlab = "Populations"))
text(x= 0.37, y = 1.5, labels = "C", cex = 2)     
text(x= 1, y= 1.36, labels= "b", cex = 2)
text(x= 2, y= 1, labels= "a", cex = 2)
text(x= 3, y= 1.15, labels= "b", cex = 2)
text(x= 4, y= 1.37, labels= "b", cex = 2)
text(x= 5, y= 0.6, labels= "a", cex = 2)
text(x= 6, y= 0.8, labels= "b", cex = 2)

title(xlab = "Populations", outer = TRUE, cex.lab= 1.7)
dev.off()

#############################################################
                 #Betadiversity calculations
#############################################################
#Data required: design; JH09_PMM_data_phyloseq_2
#Transform the count in relative abundance cpm
JH09_PMM_data_phyloseq_prop <- transform_sample_counts(JH09_PMM_data_phyloseq_2,  function(x) 1e+03 * x/sum(x))

#Cluster dendrogram
#Bray-Curtis matrix dissimilarity
BC <- phyloseq::distance(JH09_PMM_data_phyloseq_prop, "bray")
BC_cluster <- hclust(BC, "average")
plot(BC_cluster)

#Make the cluster dendrogram nicer
#First make the cluster a dendogram
dend <- as.dendrogram (BC_cluster)
plot(dend)
#Second put shapes to the populations
dend_3 <- set(dend, "leaves_pch", c(17, 18, 17, 17, 19, 18, 18, 21, 19, 2, 2,19, 19, 18, 18,2, 21, 21, 21, 17, 22, 22, 17,22,21, 21, 21, 21, 18, 18, 19, 19, 19, 22, 22, 22, 22, 22,18, 19, 17, 17, 17, 2))
plot(dend_3, ylim = c(0,0.35))
#Third remove the labels
dend_3_wolabels <-set (dend_3, "labels", c("", "", "", "", "", "", "", "","", "", "", "","", "", "", "","", "", "", "", "", "", "", "","","", "", "", "", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", ""))
plot(dend_3_wolabels, ylim = c(0,0.35))
#Fourth increase the size of the shapes
dend_3_bigger <- set(dend_3_wolabels, "leaves_cex", 1.2)

#Save the figure as a jpeg
jpeg(filename = "./Dendrogram.jpg", width = 175, height = 100, 
     units = "mm",  res=600)
par(mar=c(0.3, 4.3, 0.3, 0))
plot(dend_3_bigger, ylim = c(0,0.35), cex.axis = 1, cex.lab= 1.5, 
     ylab = "Dissimilarity")
dev.off()

#Are there significant differences between locations in bacterial community composition?
#BC distance adonis
adonis(BC ~ Location, data= design, permutations = 5000)

#############################################################
                            #ANCOM
#############################################################
#Extract OTU counts from the filtered phyoseq object
JH09_PMM_data_phyloseq_threshold  <- as.data.frame(otu_table(JH09_PMM_data_phyloseq_2))
colnames(JH09_PMM_data_phyloseq_threshold)

#Check that the order of the samples is the same in the design file and in the OTU table
design <- design[colnames(JH09_PMM_data_phyloseq_threshold), ]
rownames (design)
colnames(JH09_PMM_data_phyloseq_threshold)

#Replace the column names in the OTU table with Location ID
colnames(JH09_PMM_data_phyloseq_threshold) <- design$Location
JH09_PMM_data_phyloseq_threshold[1:10, ]

#save the table for reproducibility of the code
#write.table(JH09_PMM_data_phyloseq_threshold, file="JH09_PMM_threshold_counts_ANCOM_PC.txt", sep="\t")

# Import OTU table in 'wide format', generated manually in Excel.
JH09_PMM_threshold_wide_PC  <- read.delim("JH09_PMM_threshold_counts_ANCOM_PC_transposed.txt")

#ANCOM calculation
ancom_JH09_location_PC <- ANCOM(JH09_PMM_threshold_wide_PC, sig = 0.01, multcorr = 1)

#Generate data frame with output of the calculation
ancom_result_dataframe_PC <- as.data.frame(ancom_JH09_location_PC$detected)

#Identify the taxonomy information of the different OTUs
#Visualise ancom results with boxplot
plot_ancom(ancom_JH09_location_PC)

#Creat a list of differentially abundant OTUs
ANCOM_OTUs_PC <- as.vector(c("306710","3587343","543675","4449851","4334053","346056","New.ReferenceOTU65","246060","4394926","360440","279948","4455618","New.ReferenceOTU9","219439","New.ReferenceOTU22","513808"))

#Define the proportion of fluctuating OTUs
ANCOM_OTUs_proportion_PC <- as.data.frame (colSums (JH09_PMM_data_phyloseq_threshold[ANCOM_OTUs_PC, ]))/(colSums(JH09_PMM_data_phyloseq_threshold))*100

#write.table(ANCOM_OTUs_proportion_PC, file="ANCOM_OTUs_proportion_PC.txt", sep="\t")

#Retrieve taxonomy information of a specific OTU
dat_tax_noPlants ["306710", ]
dat_tax_noPlants ["3587343", ]
dat_tax_noPlants ["543675", ]
dat_tax_noPlants ["4449851", ]
dat_tax_noPlants ["4334053", ]
dat_tax_noPlants ["295031", ]
dat_tax_noPlants ["346056", ]
dat_tax_noPlants ["New.ReferenceOTU65", ]
dat_tax_noPlants ["246060", ]
dat_tax_noPlants ["4394926", ]
dat_tax_noPlants ["360440", ]
dat_tax_noPlants ["4455618", ]
dat_tax_noPlants ["New.ReferenceOTU9", ]
dat_tax_noPlants ["219439", ]
dat_tax_noPlants ["New.ReferenceOTU22", ]
dat_tax_noPlants ["295031", ]