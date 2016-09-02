'''
Filtering data
	  - Gather data by genome
	  - Perform standard deviation calculation
	  - Print standard deviation by genome
	  - Remove genomes below certain standard deviation

What we have
     - A dictionary named expression with values stored by cell ID and by genome
     - Want to isolate the data by line
'''

# 0. Script that extracts data from abundance.tsv files of different cell states and inputs it into a new data.tsv file
import os,sys,numpy,operator,math,scipy

from scipy import stats

def entropyCalculator(v):

    # calculating the probability distribution
    k=numpy.arange(-23-1.,10.+1.) # athan, you need to set this in the appropriate range for your values distributions. I would recommend to inspect max and min of log10 TPM, and use that to generate k. Generate like 10-20 bins from min to max values observed in your data.
    n,bins=numpy.histogram(v,bins=k)

    y=[]
    y=numpy.array(n)
    y=y/float(sum(y))

    s=scipy.stats.entropy(y)

    return s

# Now have abundace.tsv files. Next step: extract data from abudnace.tsv file. The data needing extraction is the target_id (column 1), est_counts (column 4), and tpm (column 5)
def getInfo(file, expression, cellID):
    '''
    Read data from input file
    '''
    with open(file, "r") as input_file:
        header = input_file.readline()
        for line in input_file:
            vector = line.split("\t")
            geneName = vector[0]
            tpm = float(vector[4])
            expression[cellID][geneName] = tpm

    return expression

def main():
    expression = {}
    resolution=4000
    # Input files: abundance.tsv
    inputDir='/Users/user/athan/transcriptomics/data/output/'
    cellIDs=os.listdir(inputDir) # Gathers cell ID in list
    cellIDs.sort()

    # Generates individual cell ID 
    print ('Reading input data...')
    for cellID in cellIDs:
        fileDir='%s/abundance.tsv'%(cellID) # Navigates to cell specific abundance.tsv
        path2abundance=inputDir+fileDir # Constructs abundance.tsv path
        expression[cellID] = {}
        expression = getInfo(path2abundance, expression, cellID)

    # Gather tpm data per genome
    print ('Calculating entropy...')
    analysis = {}
    sDict = {}
    maxes = []
    mins = []
    geneIDs = list(expression['71'].keys())
    geneIDs.sort()
    # Write gene IDs
    for geneID in geneIDs:
        expressionValues=[]
        # Adds tpm values for the gene cell by cell
        for cellID in cellIDs:
            value = expression[cellID][geneID]
            if value == 0:
                exp = 0
            else:
                exp = math.log(value)
            expressionValues.append(value)
        # Calculate coefficient of variance for the gene
        '''
        mean = numpy.mean(expressionValues)
        std = numpy.std(expressionValues)
        
        maximum = max(expressionValues)
        minimum = min(expressionValues)
        maxes.append(maximum)
        mins.append(minimum)
        upperRange = max(maxes)
        lowerRange = min(mins)
        print(upperRange, lowerRange)
        sys.exit()
        # Removes divide by zero error
        if sum(expressionValues) <= 25.:
            cv = 0.
        else:
            cv = std/mean
        '''
        s = entropyCalculator(expressionValues)
        sDict[geneID] = s
        
    # Sorts genes from high to low by coefficient of variance and truncates to first 5000
    variantGenes = sorted(sDict, key=sDict.get, reverse=True)[:resolution]
    for i in range(100):
        print (sDict[variantGenes[i]])

    # writing file
    print ('Writing matrix...')
    matrixFileName = '../data/Cre_Matrix_%s.txt'%str(resolution)
    with open(matrixFileName, "w") as m:
        m.write('geneName')
        for cellID in cellIDs:
            m.write('\t'+cellID)
        m.write("\n")
        for geneID in variantGenes:
            m.write('%s'%geneID)
            for cellID in cellIDs:
                value = str(expression[cellID][geneID])
                m.write('\t'+value)
            m.write('\n')
    
    return None

main()
