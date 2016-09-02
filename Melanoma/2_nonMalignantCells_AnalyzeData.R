################################################################################
### This script explores the expression of single cells from melanoma tumors ###
### and performs classification of immune cells from single cell             ###
### transcriptomes of melanoma tumors.                                       ###
################################################################################

########################
### Loading packages ###
######################## 
# Loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(scatterplot3d) # library for static 3D plotting
library(rgl) # 3D Visualization Using OpenGL
library(tictoc)
library(scatterplot3d)

##############################
### Loadata for analysis ###
##############################
# Set your working directory
setwd("/Users/user/athan/melanoma/nonMalignant/data/formatted")

# Reading the data and metadata
print('reading and treating data...')
originalData = read.csv('immuneCells.8k.genes.data.csv', header=TRUE, row.names=1)
tumorMetadata = read.csv('immuneCells.8k.genes.tumorMetadata.csv', header=TRUE, row.names=1)

# Load centroids for immune, CAF, and Endothelial cells
centroids.immune = read.csv('immuneEtcCentroids.csv',row.names=1,header=T)

# Select genes that match to centroid genes
expression = originalData[rownames(centroids.immune),]
str(expression)

##################################
### Exploring single cell data ###
##################################
# Prepare to colorize based on single cell tumor origin
colorInterpolation = colorRampPalette(brewer.pal(9,'Set1'))
col1 = colorInterpolation(length(sort(unique(tumorMetadata[,1]))))
names(col1) = sort(unique(tumorMetadata[,1]))
cols1 = as.character(col1[tumorMetadata[,1]])

##################
## Compute tSNE ##
##################

# 3D tSNE
tic()
# rtsne3d = Rtsne(t(as.matrix(expression)),dim=3,perplexity=150)
rtsne2d = Rtsne(t(as.matrix(expression)), dim=2, perplexity=150)
toc()

# Plot of 2D tSNE results
par(mar=c(4,4.5,2,1))
par(oma=c(0,0,0,0))
plot(rtsne2d$Y, main='tSNE of non-malignant cells, p = 50', col=cols1,
     pch=19, xlab = 'tSNE 1', ylab = 'tSNE2', asp = 0.7)
legend(22, 20, legend=unique(tumorMetadata[,1]), pch=19, col = col1[unique(tumorMetadata[,1])],
       bty = 'n', y.intersp = 1.2)

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols1, type='s', size=0.5, xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=unique(tumorMetadata[,1]), fill=col1[unique(tumorMetadata[,1])])
#play3d(spin3d(), duration=5)
# Need to close window by hand

################################
### Classify to immune cells ###
################################

# Classify using Correlation?
typeof(expression)
centroids.immune
cell.types <- character()
for (cellNum in colnames(expression)) {
  list1 = expression[cellNum]
  cell = as.numeric(unlist(list1))
  corList = list()
  for (cellType in colnames(centroids.immune)) {
    list2 = centroids.immune[cellType]
    type = as.numeric(unlist(list2))
    corList[cellType] <- abs(cor(cell,type))
  }
  cell.types <- c(cell.types, names(corList)[which.max(corList)])
}
cellTypeMetadata <- data.frame(cell.types)
str(cellTypeMetadata)

######################################
### Significance of predcitions by ###
### resampling of genes, and calls ###
### based on significance.         ###
###                                ###
### Cell type calls made based on: ###
###  1. p-value <= 0.05            ###
###  2. R >= 0.1                   ###
######################################
permutations = 100


###########################################
### Colorize using immune cell analysis ###
###########################################
# Prepare colors for plotting
colorInterpolation = colorRampPalette(brewer.pal(10,'Paired'))
col2 = colorInterpolation(length(sort(unique(cellTypeMetadata$cell.types))))
names(col2) = sort(unique(cellTypeMetadata$cell.types))
cols2 = as.character(col2[cellTypeMetadata$cell.types])

# 2D scatter plot
par(mar=c(4,4.5,2,1))
par(oma=c(0,0,0,0))
plot(rtsne2d$Y, main='tSNE of non-malignant cells, p = 50', col=cols2,
     pch=19, xlab = 'tSNE 1', ylab = 'tSNE2', asp = )
legend(13, 12, legend=unique(cellTypeMetadata$cell.types), pch=19, col = col2[unique(cellTypeMetadata$cell.types)],
       bty = 'n', y.intersp = 1.2)

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols2, type='s',
       size=0.5, xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=sort(unique(cellTypeMetadata$cell.types)),
         fill=col2[sort(unique(cellTypeMetadata$cell.types))])