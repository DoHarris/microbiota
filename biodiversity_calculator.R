## this caclulates the biodiversity indices of specified categories.
# the first row of the input file is for sampes (IDs).
# the second row is for the class (specified categories, e.g., CSTs), i.e., whatever trait that one desiers to access.
#Below is how the code can be used, say if one has 7 classes (e.g., CSTs) to compare 
# third row until the n-th contains the sample counts

#Differences in alpha diversity scores can be calculated outside R. For example using GraphPad

# setting working directory
setwd("name_of_your_working_directory/")

# read in data

diversity <- read.table("diversity_file.txt", header=T, sep="\t", as.is=T)

# remove cluster line
# this sign (^) functions like *
div.noc <- as.matrix(diversity[-grep("^CST",diversity[,1]),])
# change class of object to numeric for calculations
class(div.noc) <- "numeric"

# remove sample line
div.nosamp <- read.table("diversity_file.txt", header=T, skip=1, sep="\t", as.is=T)

sum.col <- colSums(div.noc)

matrix.biodiv <- matrix(0, nrow=5, ncol=ncol(div.noc), dimnames=list(c("Cluster","Simpson Index", "Dominance Index", "Shannon Index", "Shannon Equitability Index"),colnames(div.noc)))
matrix.biodiv["Cluster",] <- as.character(diversity[1,])

for(i in 1:ncol(div.noc)){
  tmp <- 0
  #tmp.eq <- 0
  tmp.shan <- 0
  for(j in 1:nrow(div.noc)){
    tmp <- tmp + div.noc[j,i]*(div.noc[j,i]-1)
    if(div.noc[j,i] != 0){
      #tmp.eq <- tmp.eq - (div.noc[j,i]/sum.col[i])*log(div.noc[j,i]/sum.col[i])
      tmp.shan <- tmp.shan + (div.noc[j,i]/sum.col[i])*log(div.noc[j,i]/sum.col[i])
    }
  }
  matrix.biodiv["Simpson Index",i] <- tmp/(sum.col[i]*(sum.col[i]-1))
  matrix.biodiv["Dominance Index",i] <- 1-(tmp/(sum.col[i]*(sum.col[i]-1)))
  matrix.biodiv["Shannon Index",i] <- -tmp.shan
  matrix.biodiv["Shannon Equitability Index",i] <- -tmp.shan/log(length(div.noc[div.noc[,i]!=0,i]))
}

## Per cluster: average + display standard deviation

# transpose matrix
biodiv.t0 <- t(matrix.biodiv)

# convert characters into numerical values
biodiv.t <- data.frame(Cluster=biodiv.t0[,1], Simpson=as.numeric(biodiv.t0[,"Simpson Index"]), 
                       Dominance=as.numeric(biodiv.t0[,"Dominance Index"]), 
                       Shannon=as.numeric(biodiv.t0[,"Shannon Index"]), 
                       Shannon_Equitability=as.numeric(biodiv.t0[,"Shannon Equitability Index"]))
                       
# aggregate per cluster to mean and standard deviation
test <- aggregate(apply(biodiv.t[,-1],2,as.numeric), # convert to numeric data
                  by=list(biodiv.t[,1]), # collapse per cluster name
                  function(x)paste("median:",round(median(x),4),", mean:",round(mean(x),4),", sd:", round(sd(x),4), sep="")) # functions that calculates average and standard deviation and displays them with text



## boxplot

library(ggplot2)
library(reshape2)

# preparing table

# change rownames to sample + cluster
#rownames(biodiv.t) <- paste(gsub("X","",rownames(biodiv.t)), biodiv.t[,"Cluster"], sep="_")
biodiv.t$samples_clusters <- paste(gsub("X","",rownames(biodiv.t)), biodiv.t[,"Cluster"], sep="_")


# remove Cluster column
#biodiv.mat <- biodiv.t[,-1]

# convert to the "long" format
biodiv.long <- melt(biodiv.t)

# split sample name and clusters
biodiv.df <- data.frame(biodiv.long, 
                        samples=unlist(lapply(strsplit(as.vector(biodiv.long$samples_clusters), split="_"),function(x)x[[1]])),
                        clusters=unlist(lapply(strsplit(as.vector(biodiv.long$samples_clusters), split="_"),function(x)x[[2]])))

# name remaining columns
colnames(biodiv.df)[3] <- "indices"

# convert value column into numeric 
class(biodiv.df$value) <- "numeric"

# boxplot for means

# with outliers, add "outlier.shape=NA" 
# fill=clusters set to color the boxes
# http://stackoverflow.com/questions/10805643/ggplot2-add-color-to-boxplot-continuous-value-supplied-to-discrete-scale-er
# http://docs.ggplot2.org/current/geom_boxplot.html
# http://htmlcolorcodes.com/color-chart/web-safe-color-chart/

p <- ggplot(biodiv.df, aes(x=indices, y=value)) + geom_boxplot(lwd=1.5, aes(color=clusters), position=position_dodge(1)) + 
  scale_colour_manual(labels = c("CST-1", "CST-2", "CST-3", "CST-4", "CST-5", "CST-6", "CST-7", "CST-8"), 
                      values = c("#0000FF", "#00FF00", "#00FFFF", "#006600", "#CC00CC", "#FFFF33", "#FF0000"))

# add vertical lines
p1 <- p + geom_vline(xintercept=c(1.5, 2.5, 3.5), linetype="dashed", color="black", size=0.5)



#http://www.cookbook-r.com/Graphs/Legends_%28ggplot2%29/

p1 + theme(legend.title = element_text(colour="black", size=20, face="bold"), # Title appearance
          legend.text = element_text(colour="black", size = 18), text=element_text(size=28)) # Label appearance


#This written by Sarah Bonnin and Harris Onywera
#June 2016
