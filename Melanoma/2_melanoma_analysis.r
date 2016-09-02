#################################################################
# @Course: Systems Biology of Disease                           #
# @Rscript: dimensionalityReduction.R                           #
# @Author: Adrian de Lomana and Chris Plaisier                  #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

###
### This script performs data exploration and classification of melanoma single cell transcriptomes.
### Specifically it works with transcriptomes of malignant and non-malignant cells.
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(devtools)
library(ggbiplot)
library(caret)
library(e1071)
library(ggfortify)
library("Rtsne")
library(lattice)
library('scatterplot3d')
library(RColorBrewer)
library(rgl)

# 0.2. defining working directory
setwd('/Users/user/athan/melanoma/malignantCells/data/formatted')

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
mydata <- read.csv("malignant.8k.genes.data.csv", header = TRUE, row.names = 1)
mymetadata = read.csv("malignant.8k.genes.tumorMetadata.csv", header = TRUE, row.names = 1)
#str(mydata)
#str(mymetadata)
td = t(mydata)
inputDataMatrix = as.matrix(td)

# 1.1 Tagging cell names with tumor labels
cellnames = list()
x<-1
for (cell.name in colnames(mydata)) {
  cellnames[[x]] <- cell.name
  tumorlabel <- mymetadata[cell.name,]
  attr(cellnames[[x]], "Tumor") <- tumorlabel
  x <- x + 1
}
str(cellnames)


# Returns cancer label
tumortypes <- levels(attributes(cellnames[[8]])$Tumor)

# Setting plot colors
i <- 1
for (cellname in cellnames) {
  if (attributes(cellnames[[i]])$Tumor == 'Mel78')  {
    attr(cellnames[[i]], "Color") <- 'red'  
  }
  else if (attributes(cellnames[[i]])$Tumor == 'Mel79')  {
    attr(cellnames[[i]], "Color") <- 'orange'  
  }
  else if (attributes(cellnames[[i]])$Tumor == 'Mel80')  {
    attr(cellnames[[i]], "Color") <- 'yellow'  
  }
  else if (attributes(cellnames[[i]])$Tumor == 'Mel81')  {
    attr(cellnames[[i]], "Color") <- 'green'  
  }
  else if (attributes(cellnames[[i]])$Tumor == 'Mel84')  {
    attr(cellnames[[i]], "Color") <- 'blue'  
  }
  else if (attributes(cellnames[[i]])$Tumor == 'Mel88')  {
    attr(cellnames[[i]], "Color") <- 'purple'  
  }
  i <- i + 1
}

# 2. dimensionality reduction analysis
# 2.0. setting some variables for plotting
tumorLabels=as.character(mymetadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)

# 2.1. PCA
                                        
# 2.2. t-SNE
# 2.2.1 2d low resolution of tSNE
tsneResult2d = Rtsne(inputDataMatrix, dim=2, theta=0.5, verbose=TRUE, perplexity=5)
plot(tsneResult2d$Y, main='tSNE of malignant cells, p = 5', col=plottingColors[tumorLabels],
     pch=19, xlab = 'tSNE 1', ylab = 'tSNE2')
par(mar=c(4,4.5,2,1))
par(oma=c(0,0,0,0))
legend(47, 47, legend=unique(tumorLabels), pch=19, col = plottingColors[unique(tumorLabels)],
       bty = 'n', y.intersp = 1.2)

#2.2.2 2d high resolution of tSNE
tsneResult2d = Rtsne(inputDataMatrix, dim=2, theta=0.5, verbose=TRUE, perplexity=50)
tumorLabels=as.character(mymetadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
plot(tsneResult2d$Y, main='tSNE of malignant cells, p = 50', col=plottingColors[tumorLabels],
     pch=19, xlab = 'tSNE 1', ylab = 'tSNE2')
par(mar=c(4,4.5,2,1))
par(oma=c(0,0,0,0))
legend(18, -11, legend=unique(tumorLabels), pch=19, col = plottingColors[unique(tumorLabels)],
       bty = 'n', y.intersp = 1.2)



# 2.2.2. 3d high resolution of t-SNE
tsneResult3d = Rtsne(inputDataMatrix,dim=3,theta=0.5,verbose=TRUE,perplexity=50)
par3d(windowRect=c(50,50,700,700))
plot3d(tsneResult3d$Y, col=plottingColors[tumorLabels], type='s', size=0.5, xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=unique(tumorLabels), pch=19, col = plottingColors[unique(tumorLabels)],
         bty = 'n', y.intersp = 1.2)
# ecb = function(x,y){ 
#   plot(x, t = 'n', main = 'tSNE of malignant melanoma cells',
#                           xlab = 'tSNE 1',
#                           ylab = 'tSNE 2');
#   text(x, labels=tumortypes, col = colors[tumortypes])
#   }
# tsne_mel = tsne(tsnematrix[1:5,], epoch_callback = ecb, perplexity=50)

# plotting

                                        
# 2.2.2. high resolution of t-SNE

# 2.2.3. high resolution of t-SNE showing 3D

                                        
# interactive 3D scatter plot. You will need to close the interactive plot window when you're done with it



