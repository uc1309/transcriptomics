# 0. Script that extracts data from abundance.tsv files of different cell states and inputs it into a new data.tsv file
import os,sys

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

def simplify():
    '''
    Removes extraneous data from output file
    '''

def main():
    expression = {}
    analysis = {}
    # Input files: abundance.tsv
    inputDir='/Users/user/athan/transcriptomics/data/output/'
    cellIDs=os.listdir(inputDir) # Gathers cell ID in list                        

    # Generates individual cell ID 
    print ('Reading...')
    for cellID in cellIDs:
        fileDir='%s/abundance.tsv'%(cellID) # Navigates to cell specific abundance.tsv
        path2abundance=inputDir+fileDir # Constructs abundance.tsv path
        expression[cellID] = {}
        expression = getInfo(path2abundance, expression, cellID)

    print ('Writing...')
    outputFileName = 'output.txt'
    geneIDs = expression[cellID].keys()
    # Write header with cell IDs
    with open(outputFileName, "w") as output_file:
        output_file.write("gene_ID")
        for cellID in cellIDs:
            output_file.write(cellID + "\t")
        output_file.write("\n")
        # Write gene IDs
        for geneID in geneIDs:
            output_file.write(geneID)
            # Write tpm value
            for cellID in cellIDs:
                value = str(expression[cellID][geneID])
                output_file.write(value + "\t")
            output_file.write('\n')

    return None

main()




