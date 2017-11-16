# Thanks for reading this!

# This is a collection of my R codes that I think should be helpful for your research.
# In general, what is below are functions that you can load into your own script.
# These functions include from loading the data to perform simple analysis.
# They are designed for amplicon sequencing data.

# I basically want to put all of them down here is to shorten my thinking time
# from a "need" to a "solution",
# so consider this as a reference book. Of course, you don't need to follow everything.
# Eventually, you will develop you own way of coding (and not just in R!).

# Zewei Song (songzewei@genomics.cn)

# You can summon this script in your code by includig the following lines:
# function_path <- 'C:/SongLab/Tools/R_Scripts/Amplicon_Data' # This is the path of this script in your computer
# setwd(function_path)
# source('ampliconR.R')

# A better way (for me) is to keep this one open in RStudio, you can select all (or the one you need)
# functions and Ctrl+Enter to run them.

# The example data is an OTU table, and a treatment table
# The OTU table contains some samples that are missing in the treatment table,
# While the treatment table contains some samples that are missing in the OTU table, and
# are in a different order.
# The data is a subset of a larger table, so the abundance of some OTUs are all zero.
# What a mess! And we are going to deal with this.

# You can now select all the functions and run them (this load them into memory),
# all the way down to the long section line (##### .... ###).
# Then just follow the code in the example section to the end.
# In RStudio, functions and sections can be closed, open them use the triangle by the line number.

### FUNCTIONS #######################################################
# Read in the OTU table (with samples in columns and OTUs in rows) from a tab delimited file as default
# Use meta_column to set the number of meta data columns. You need to make sure the number is right.
# The function use this number to decide where to trim the OTU table.
# This number default to 1 so you have one additional column for taxa.
read_otu_table <- function(otu, meta_column = 1, delimiter = '\t', remove_zero = FALSE){
  print(paste0('Reading in ', otu, ' ...'))
  otu_table <- read.csv(otu, header=TRUE, row.names=1, check.names = FALSE, sep = delimiter)
  print(dim(otu_table)[2])
  otu_table <- data.frame(t(otu_table[,1:(dim(otu_table)[2]-meta_column)]))
  sample_number = dim(otu_table)[1]
  otu_number = dim(otu_table)[2]
  print(paste0('Found ', sample_number, ' samples and ', otu_number, ' OTUs.'))
  print(paste0('This OTU table contains ', count_zero(otu_table), ' OTUs with zero abundance.'))
  if(remove_zero){
    #print('Zero abundance OTUs were removed from the data frame.')
    return(remove_zero(otu_table))
  }
  else{
    print('Zero abundance OTUs were kept in the data frame.')
    return(otu_table)
  }
}

# Count the number of rows that all elements are zero
# You cannot have non-abundance column in the OTU table (i.e. need to make sure meta_column in read_otu_table() is right)
count_zero <- function(data){
  a <- data.frame(t(data))
  a$rowSum <- data.frame(rowSums(a))
  a <- a[(a$rowSum == 0),]
  return(dim(a)[1])
}

# Remove rows if all elements are zero from a data frame
remove_zero <- function(data){
  #print(paste0('This OTU table contains ', dim(data)[2], ' OTUs.'))
  temp_row <- row.names(data) # Save the right row names (i.e. sample names)
  a <- data.frame(t(data)) # Data is transposed
  a$rowSum <- data.frame(rowSums(a))
  a <- a[order(a$rowSum, decreasing=TRUE),]
  a <- a[(a$rowSum>0),]
  a <- subset(a, select = -(length(colnames(a))))
  a <- data.frame(t(a)) # Data is transposed back
  row.names(a) <- temp_row # R replace all "-" to ".", this fix it.
  print(paste0('The current OTU table has ', dim(a)[2], ' none-zero OTUs, they were returned to the new table.'))
  return(a)
}

# Read in experimental treatment
read_treatment <- function(treatment, delimiter = '\t'){
  print(paste0('Reading in ', treatment, ' ...'))
  trt <- read.csv(treatment, header=TRUE, row.names=1, check.names=FALSE, sep=delimiter)
  print(paste0('Found ', dim(trt)[1], ' samples with ', dim(trt)[2], ' treatment.'))
  print(str(trt))
  return(trt)
}

# Grab treatment meta data using the sample names in the OTU table
# This can result in a subset of treatment that is read in by read_treatment()
# This will also fix the problem if the order of the sample in treatment is not the same
# as the one in the OTU table.
get_treatment <- function(otu, trt){
  for (sample in row.names(otu)) {
    if (! sample %in% row.names(trt)){
      print(paste0('Sample ', sample, ' is not in the treatment data, please check.'))
    }
  }
  i <- intersect(row.names(otu), row.names(trt))
  return(trt[i,])
}

# Read in the dataset, including OTU table and treatment
# This will create a data list with otu and trt, probabaly can add dist to it in the future.
get_data <- function(otu, trt, remove_zero=TRUE){
  trt <- get_treatment(otu, trt) # get the intersect between OTU and Treatment
  otu <- remove_zero(otu[row.names(trt),])
  return(list(otu = otu, trt = trt))
}

# Calculate the distance/dissimilarity matrix
add_dist <- function(data_list, transform='total', dist='bray'){
  d <- vegdist(decostand(data_list$otu, method=transform), method=dist)
  data_list[['dist']] <- d
  return(data_list)
}

# Convert a dist Class to column based format
dist2column <- function(dist){
  b <- data.frame(t(combn(rownames(as.matrix(dist)),2)), as.numeric(dist))
  names(b) <- c('c1', 'c2', 'distance')
  return(b)
}

# Convert the column based distance matrix back to a dist Class
column2dist <- function(column){
  library(reshape)
  col_names <- colnames(column)
  samples <- unique(c(as.vector(column[,col_names[1]]), as.vector(column[,col_names[2]])))
  
  c <- reshape(column, direction='wide', v.names=col_names[3], timevar=col_names[2], idvar=col_names[1])
  c <- data.frame(c, row.names = 1)
  
  # Now we need to add the missing column and row for the distance matrix
  # Create empyty vectors
  new_col <- rep(NA, dim(c)[1]) # New column vector
  new_col[1] <- 1.0
  c <- cbind(new_col, c)
  new_row <- rep(NA, dim(c)[2]) # New row vector
  new_row[length(new_row)] <- 1.0
  c <- rbind(c, new_row)
  row.names(c) <- samples # Add the right sample names
  colnames(c) <- samples
  c[lower.tri(c)] <- t(c)[lower.tri(c)] # Cope the upper triangle values to lower triangle
  d <- as.dist(c) # Convert the matrix the dist Class
  
  return(d)
}

# Give a dist object, slice based on the provided vector (need to have the same sample name)
slice_dist <- function(d, s){
  a <- as.data.frame(as.matrix(d)) # convert the dist to data frame
  a <- a[s,s] # Slice the data frame base on the Bolean vector
  b <- as.dist(a) # Convert back to dist object
  return(b)
}

### EXAMPLES ##################################################################
# Set your working directory
rm(list = ls()); # Use this if you have too many stuff in your memory/environment, but don't forget to load the functions first
function_path <- 'C:/SongLab/Tools/R_Scripts/Amplicon_Data' 
setwd(function_path)
library(vegan) # this is the community ecology package we need. Take a guess why it is named "vegan" :D

# Read in the OTU table
otu_file = 'amplicon_example_otu.txt'
otu <- read_otu_table(otu_file, meta_column = 1) # This is a tab-delimited file with one taxonomic column
# You can see that 6217 out of 7222 are all-zeros, but we kept all of them at this point.

# Read in the treatment
trt_file = 'amplicon_example_treatment.txt'
trt <- read_treatment(trt_file)
# We have 12 samples in the treatment data, but only 11 samples in the OTU table.

# Combine the OTU table and treatment into a single data list
# And removed all blank OTUs
# I use dat instead of data since data is a popular parameters for many functions. You can use data, but it is easier to create confusion.
dat <- get_data(otu, trt, remove_zero = TRUE)

# I found it more convenient to put all data into a single list,
# since some times we need to work with several OTU tables.
# We can now add a Bray-Curtis dissimilarity matrix (or any kind) into the list
# For fungal data, I always do a Hellinger trasnform to deal with the highly skew distribution.
dat <- add_dist(dat, transform = 'hellinger', dist = 'bray')

# Let's see how the dist matrix fold out on a NMDS plot with two dimensions
m <- metaMDS(dat$dist, method = 'bray', trymax = 100) # We will get a warning sinc we only have 9 samples
stressplot(m)
ordiplot(m, type="text", display = "site", main="Bray-Curtis")
score = data.frame(scores(m, display = 'sites')) # This is the coordination that you can use to plot
write.csv(score, file='nmds_score.csv') # I prefer to write them to file, and draw plot in other softwares.

# Analyze a subset of your data
# The dist matrix can be futher sliced into subset.
# I prefer to slice the dist matrix instead of original data so you don't need to calculate
# the matrix again. This also makes your code simpler.
pcr1 <- dat$trt$PCR == 'pcr1' # We use the treatment to create a slicer for all samples coming from PCR1
# You can combine the logic operator (!, &, |, ...) for complicated slicers.
# pcr1 has the same length as the number of samples in our data, with either TRUE or FALSE.
# Let's do another NMDS with the subset pcr1
m_pcr1 <- metaMDS(slice_dist(dat$dist, pcr1), method = 'bray')
stressplot(m_pcr1)
ordiplot(m_pcr1, type="text", display = "site", main="Bray-Curtis")
# If you want, you can also write a small function to wrap all these NMDS steps together.

# Multivariate analysis
# We want to evaluate the effect size of Plot and PCR, and their interactions.
adonis(dat$dist ~ Plot * PCR, data = data$trt)
# For real (larger) OTU tables, you will get P values, but not for this example.
# And of course we can use slicer here, too. But now we need to slice the treatment data too.
adonis(slice_dist(dat$dist, pcr1) ~ Plot, data = dat$trt[pcr1,])

# If you want to explore, here are more functions that you can try by yourself.
# Take a look at http://mb3is.megx.net/gustame if you are not sure what to do
require(ape)
?pcoa #PCoA

?envfit # Fit environmental factors
?anosim # similar to ADONIS, but only compare among groups
?cca # A constrained multivariate analysis
?dbrda # distance based redundancy analysis

### STOP HERE ##################################################################
# Code below is something I'm still testing. Feel free to read them. Don't panic.
# Example for dist2column
# Given a distance matrix (a) in Class 'dist' (as generated by vegan:vegdist), generated from the OTU table x
# Use this command to convert it to column based format (b)
working_path = 'E:/GoogleDrive/Research/Tools/scripts'
setwd(working_path)
example_otu_table = 'R_example_otu_table.txt'
otu_table <- read.csv(example_otu_table, header=TRUE, row.names=1, check.names = FALSE, sep = "\t")
otu_table <- t(otu_table)

library(vegan)
a <- vegdist(otu_table) # calculate the distance matrix (or dissimilarity matrix)
b <- data.frame(t(combn(rownames(as.matrix(a)),2)), as.numeric(a)) # convert a to column based format
names(b) <- c('c1', 'c2', 'distance') # Add column names to b

# Example for column2dist
#I can use reshape or reshape2 to convert the column based format back to matrix
# The first two columns need to be c1 and c2
library(reshape2)
# get the sample names from the column data b:
samples <- unique(c(as.vector(b$c1), as.vector(b$c2)))
col_names <- colnames(b)
samples <- unique(c(as.vector(b[,col_names[1]]), as.vector(b[,col_names[2]])))
c <- reshape(b, direction='wide', v.names='distance', timevar='c2', idvar='c1')
c <- data.frame(c, row.names = 1)
# Now we need to add the missing column and row for the distance matrix
# Create empyty vectors
new_col <- rep(NA, dim(c)[1]) # New column vector
new_col[1] <- 1.0
c <- cbind(new_col, c)
new_row <- rep(NA, dim(c)[2]) # New row vector
new_row[length(new_row)] <- 1.0
c <- rbind(c, new_row)
row.names(c) <- samples # Add the right sample names
colnames(c) <- samples
c[lower.tri(c)] <- t(c)[lower.tri(c)] # Cope the upper triangle values to lower triangle
d <- as.dist(c) # Convert the matrix the dist Class