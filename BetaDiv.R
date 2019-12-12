#Read in data
#set working directory - where you keep your biom and metadata file
#-----------------------------------------------------------------------------------------------------
#This script requires a metadata file (e.g., mapp.txt), with the first column being "Sample_IDs" and the subsequent columns being the other particulars about the participant.
#In our analysis, the particulars included the following: HRHPV, HIV, HIVHRHPV, CD4, HPV, HIVHPV, CST, and CST_pooled.

#The other file required is a rarefied OTU table, e.g, otu_table_even13014.tax.biom

setwd("name_of_working_directory/")

#Load in the phyloseq & other packages
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("Heatplus")
biocLite("Rcpp")
biocLite("biom")


###One might get an error when installing the Bioconductor packages, depending on the R version that he/she is using.
###If so, install Bioconductor packages using BiocManager; see https://bioconductor.org/install

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")




install.packages("phyloseq") #In case you find the error that says package 'phyloseq' is not available, use the folowing command BiocManager::install():
BiocManager::install("phyloseq")


install.packages("ape")
install.packages("BAT")
install.packages("codetools")
install.packages("cluster")
install.packages("corrplot")#for cor.mtest
install.packages("dplyr")
install.packages("fifer")
install.packages("ggplot2")
install.packages("graphics")
install.packages("matrixStats")#rowSds
install.packages("NMF")#for heatmap function
install.packages("metagenomeSeq")#differential abundance testing
install.packages("pheatmap")
install.packages("plyr")
install.packages("psych")#corr.test
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("scales")
install.packages("vegan")
install.packages("picante")
install.packages("gridExtra")
install.packages("knitr")
install.packages("colorspace")
install.packages("Biostrings")
install.packages("permute")
install.packages("lattice")
install.packages("nlme")
install.packages("gplots")
install.packages("Heatplus")
install.packages("RColorBrewer")
install.packages("biom")
install.packages("devtools")




library(phyloseq)
library(ape)
library(BAT)
library(cluster)
library(corrplot)#for cor.mtest
library(dplyr)
library(fifer)
library(ggplot2)
library(graphics)
library(matrixStats)#rowSds
library(NMF)#for heatmap function
library(metagenomeSeq)#differential abundance testing
library(pheatmap)
library(plyr)
library(psych)#corr.test
library(reshape2)
library(RColorBrewer)
library(scales)
library(vegan)
library(picante)
library(gridExtra)
library(knitr)
library(colorspace)
library(Biostrings)
library(permute)
library(lattice)
library(nlme)
library(gplots)
library(Heatplus)
library(RColorBrewer)
library(biom)
library(devtools)


#IMPORT DATA (OTU table generated in QIIME, now in the biom format)
import_biom(BIOMfilename = "otu_table_even13014.tax.biom", verbose = TRUE)
phy <- import_biom(BIOMfilename = "otu_table_even13014.tax.biom", verbose = TRUE)

# Number of OTUs in the phyloseq object
ntaxa(phy)

#IMPORT METADATA
meta <- read.table("mapp.txt", sep = "\t", header =TRUE, row.names=1)
head(meta)

#Check if the same sample numbers in meta file as the phy object
head(sample_names(phy))
length(sample_names(phy))#198
length(rownames(meta))

## Add metadata to phyloseq object.
##assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta
##Check sample numbers in the merged phy object.
nsamples(phy)

#phy function features
phy
otu_table(phy)
sample_data(phy)
tax_table(phy)
sample_names(phy)
sample_variables(phy)
taxa_names(phy)
nsamples(phy)

sample_names(phy)
sample_data(phy)$HRHPV
sample_data(phy)$HIV
sample_data(phy)$HIVHRHPV
sample_data(phy)$CD4
sample_data(phy)$HPV
sample_data(phy)$HIVHPV
sample_data(phy)$CST
sample_data(phy)$CST_pooled


#Data clean up
colnames(tax_table(phy))
#Replace "Rank1" with Kingdom etc,
colnames(tax_table(phy)) <-  c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus", "Species")
# Clean taxonomic annotations
tax_table(phy)[,"Kingdom"] <- sub("k__","",tax_table(phy)[,"Kingdom"])
tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])


#(check that the sample names match in all cases)
length(intersect(rownames(meta),sample_names(phy)))


#MORE phy function features
meta
sample_data(phy)$HRHPV
sample_data(phy)$HIV
sample_data(phy)$HIVHRHPV
sample_data(phy)$CD4
sample_data(phy)$HPV
sample_data(phy)$HIVHPV
sample_data(phy)$CST
sample_data(phy)$CST_pooled

#AGAIN (check that the sample names match in all cases)
rownames(meta)
length(rownames(meta))
sample_names(phy)
length(sample_names(phy))
length(intersect(rownames(meta),sample_names(phy)))

#MERGING  OTU SCRIPTS
#joey711.github otu merging
  # How many taxa before agglomeration?
  ntaxa(phy)

  # agglomerate at the Genus taxonomic rank
  phy <- tax_glom(phy, taxrank="OTU") 

  # How many taxa after agglomeration?
  ntaxa(phy)
  
  # print the available taxonomic ranks. Shows only 1 rank available, not useful for tax_glom
  colnames(tax_table(phy))
  
#ADD METADATA TO PHYLOSEQ OBJECT
#-------------------------------------------
sample_data(phy) <- meta
nsamples(phy)#144
  
H <-phy

#EXPLORATORY....
H <- subset_taxa(H, taxa_sums(H)>0)#keep only taxa where there are positive values (i.e. not only zeros)
ntaxa(H)# 5859

#EXPLORE NUMBER OF READS PER SAMPLE AND DISTRIBUTIONS:
readsumsdf = data.frame(ncounts = sort(taxa_sums(H), TRUE), sorted = 1:ntaxa(H), 
                        type = "otu")
readsumsdf = rbind(readsumsdf, data.frame(ncounts = sort(sample_sums(H), TRUE), 
                                          sorted = 1:nsamples(H), type = "Samples"))
readsumsdf

pdf("RPlots/Read_counts.pdf")
title = "Total number of counts"
p = ggplot(readsumsdf, aes(x = sorted, y = ncounts)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()

ntaxa(H)
nsamples(H)
min(sample_sums(H))
max(sample_sums(H))
H


#generate phylogenetic tree
random_tree = rtree(ntaxa(H), rooted=TRUE, tip.label=taxa_names(H))
H.tree = merge_phyloseq(H,random_tree)
H.tree
nsamples(H.tree)


####In case of alpha diversity analyses, one do as folows:
####PLOT RICHNESS AFTER NORMALISATION
###ALPHA DIVERSITY
##SHANNON
#define ggplots: categories (e.g., CD4) = <350 vs. >350

str(H.tree)
nsamples(H.tree)

###Plots is a pre-existing directory
pdf("Plots/Shannon_Observed_Chao1.pdf")
p1 <- plot_richness(H.tree,x = "CD4",color = "CD4",measures=c("Shannon","Observed","Chao1"), 
                    title = paste0("Standardized to total reads, N=",nsamples(H.tree)))
p1
p1 + geom_boxplot()
dev.off()

pdf("RPlots/Shannon.pdf")
p2 <- plot_richness(H.tree,x = "CD4",color = "CD4",measures=c("Shannon"), 
                    title = paste0("Standardized to total reads, N=",nsamples(H.tree)))
p2
p2 + geom_boxplot()
dev.off()


#For beta diversity analyses, one does as folows:
##BETA DIVERSITY
#####RPLOT ORDINATION
#MULTIDIMENSIONAL SCALING
#BRAY
pdf("RPlots/Bray_total.pdf")
GP.ord <- ordinate(H.tree, "MDS", "bray")
BR.1.1 = plot_ordination(H.tree, GP.ord, type = "samples", color = "CD4", shape = "HIV")
BR.1.1  = BR.1.1 + geom_point(size = 4) + ggtitle("BRAY: MDS of rarefied OTUs")+theme(axis.text=element_text(size=14, face="bold"),
                                                                           axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14))+labs(shape="HPV", colour="HIV")
BR.1.1
dev.off()

#MDS ("PCoA") on UNIFRAC DISTANCES
#phylogenetic tree
H.tree

##UNIFRAC
#https://www.rapidtables.com/web/color/RGB_Color.html

#These colours were used for the HIV, CD4, HRHPV, and HPV.
Harris_palette <- c("HPV+"="#00FFEF","HPV-"="#00FF00")

#These colours were used for HIVHRHPV and HIVHPV.
Harris_palette <- c("HIV+HPV+"="#00FFEF","HIV+HPV-"="#FFBF00","HIV-HPV+"="#006400","HIV-HPV-"="#00FF00")

#These colours were used for CST_pooled.
Harris_palette <- c("CSTs 2-6"="#000000","CST-1"="#00FF00")

#These colours were used for CST.
Harris_palette <- c("CST-1"="#00FF00","CST-2"="#FFF200","CST-3"="#0000FF","CST-4"="#FFBF00","CST-5"="#FF0000","CST-6"="#FF00FF")




par(pty="s")
pdf("RPlots/W_Unifrac_AAA.pdf")
ordu_1a = ordinate(H.tree, "PCoA", "unifrac", weighted=TRUE)
ordu_1b = ordinate(H.tree, "PCoA", "unifrac", weighted=FALSE)
PCOA1a = plot_ordination(H.tree, ordu_1a, type = "samples", color = "CST", title = "W-Unifrac: PCoA of OTUs")
PCOA1a + geom_point(size = 5)  +  scale_colour_manual(values = Harris_palette) + coord_fixed(ratio = 1) + theme(text = element_text(size=20))


dev.off()


pdf("RPlots/UW_Unifrac_CD4.pdf")
PCOA1b = plot_ordination(H.tree, ordu_1b, type = "samples", color = "CD4", shape = "CD4", title = "UW-Unifrac: PCoA of OTUs")
PCOA1b + geom_point(size = 4) +  scale_colour_manual(values = Harris_palette)
dev.off()


#Differences in beta diversity measures (between sample categories, e.g., CD4<350 vs. CD4>350) were tested using multivariate analysis.

#Permanova1 
H_UW <- phyloseq::distance(H.tree, method = "unifrac", weighted=FALSE)
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(H.tree))
# Adonis test
adonis(H_UW ~ CD4, data = sampledf)



#Permanova2 
H_W <- phyloseq::distance(H.tree, method = "unifrac", weighted=TRUE)
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(H.tree))
# Adonis test
adonis(H_W ~ CD4, data = sampledf)


#This written by Jerome Wendoh and Harris Onywera
#December 2016
