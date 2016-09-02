import sys
import numpy as np
import matplotlib
import matplotlib.pyplot
import matplotlib.patches as mpatches
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# Load matrices for data and labels
resolution=2150
print ('Loading data matrix...')
# X = np.loadtxt("/Users/user/athan/transcriptomics/data/Cre_Matrix_%s.txt"%str(resolution), skiprows=1,usecols=range(1,26))
X = np.loadtxt("/Users/user/Downloads/selectedMatrix.2fold.combined.1355.2fold.plus.794.30minsTransient.txt", skiprows=1, usecols=range(1,16))
X_transposed = np.transpose(X)
print(X_transposed)

# associating labels to colors

plotInfo = {}
theColors =['b','r','m', 'y', 'g', 'c', 'k', 'w']
theMarkers=['o','s','^','v']
print ('Loading data labels...')


# Generates TSNE model
print ('Generating tSNE model...')
model = TSNE(n_components=2)

# Generates PCA models
print ('Generating PCA model...')
#pca = PCA(n_components=2, copy=True)

# computes t-SNE 
print ('Computing on data...')
Z = model.fit_transform(X_transposed)

# Computes PCA
# Z = pca.fit_transform(X_transposed)
#print(pca.explained_variance_ratio_)

# Plot visualization
print ('Plotting data...')

'''
reds = [0, 1, 11, 12, 19, 'Time: 0']
yellows = [5, 'Time: 12m'] # 12 m
greens = [10, 'Time: 60m'] # 60 m
cyans =[8, 13, 14, 'Time: 30m'] # 30 m
blues = [15, 16, 'Time: 4h'] # 4 h
magentas = [17, 18, 'Time: 8h'] # 8 h
'''

reds = [0, 1]
yellows = [2, 3]
greens = [4, 5]
cyans = [6, 7]
blues = [8, 9]
magentas = [10, 11]
blacks = [12, 13]

'''
firstLongTimeSingles = [12, 14, 16, 18]
firstLongTimePairs = [11, 13, 15, 17]
secondLongTimes = [19, 20, 21, 22, 23, 24]
shortTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
'''

shortTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
legends = [0, 2, 4, 6, 8, 10, 12, 14]

sampleIDs = range(len(Z))

for sampleID in sampleIDs:
    print (sampleID)
    if sampleID in reds:
        color = theColors[1]
        labelTime = 'Time: 2-4m'
    elif sampleID in magentas:
        color = theColors[2]
        labelTime = 'Time: 4-8h'
    elif sampleID in yellows:
        color = theColors[3]
        labelTime = 'Time: 8-12m'
    elif sampleID in greens:
        color = theColors[4]
        labelTime = 'Time: 18-24m'
    elif sampleID in cyans:
        color = theColors[5]
        labelTime = 'Time: 30-45m'
    elif sampleID in blues:
        color = theColors[0]
        labelTime = 'Time: 1-2h'
    else:
        color = theColors[6]
        labelTime = 'Time: 12-48h'

    '''
    if sampleID in firstLongTimeSingles:
        marker = theMarkers[0]
        labelRun = 'First Long Time Run'
    elif sampleID in firstLongTimePairs:
        marker = theMarkers[1]
        labelRun = 'First Long Time Run'
    elif sampleID in secondLongTimes:
        marker = theMarkers[2]
        labelRun = 'Second Long Time Run'
    elif sampleID in shortTimes:
        marker = theMarkers[3]
        labelRun = 'Short Time Run'
    '''

    if sampleID in shortTimes:
        marker = theMarkers[0]
    
    print (color, marker)
    if sampleID in legends:
        matplotlib.pyplot.plot(Z[sampleID][0],Z[sampleID][1],'o', color=color, marker=marker, label='%s'%(labelTime))
    else:
        matplotlib.pyplot.plot(Z[sampleID][0],Z[sampleID][1],'o', color=color, marker=marker)
    

'''    
for i in range(len(Y)):
    matplotlib.pyplot.plot(Y[i][0],Y[i][1], 'ok')
'''

'''
red_patch = mpatches.Patch(color='red', label='Time: 0')
matplotlib.pyplot.legend(handles=[red_patch])
matplotlib.pyplot.show()
'''
matplotlib.pyplot.legend(loc=(0.75, 0.45))
matplotlib.pyplot.ylabel('tSNE 2')
matplotlib.pyplot.xlabel('tSNE 1')
matplotlib.pyplot.title('tSNE of N-starved C. reinhardtii (2150 genes)')
#matplotlib.pyplot.xlim(-80,120)
#matplotlib.pyplot.ylim(-30,60)

matplotlib.pyplot.savefig("tSNEResult_%s.pdf"%str(resolution))
#matplotlib.pyplot.savefig("tSNEResult_%s.png"%str(resolution))

