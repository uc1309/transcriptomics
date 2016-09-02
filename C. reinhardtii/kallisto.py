# Script to analyze reads using kallisto
import os,sys

# User defined variables
fastqDir='/Users/user/athan/transcriptomics/data/cleanfastq/'
genomeIndex='/Users/user/athan/transcriptomics/data/genomeIndex/Creinhardtii_281_v5.5.idx'
kallistoDir='/Users/user/athan/transcriptomics/data/output/'

# User defined methods

# Identify important file data
def fileID(file):
    fileNum=file.split('_')[0][7:] # File number
    fileDir=file.split('_')[1][:] # File direction
    fileType=file.split('.')[0][12:] # File type (paired or unpaired)
    fileTag=[fileNum, fileDir, fileType]
    return fileTag
    
# Select input files
fastqFiles=os.listdir(fastqDir)

# Define single end samples
singleEndSamples=['81', '83', '85', '87']

# Execute kallisto command
for fastqFile in fastqFiles:
    path2fastq=fastqDir+fastqFile

    lastNumbers=fastqFile.split('.')[0][-2:]
    if lastNumbers in singleEndSamples:
        command='kallisto quant -i %s -o %s%s --single -l 200 -s 20 %s'%(genomeIndex,kallistoDir,lastNumbers,path2fastq)
        os.system(command)
    else:
        # Sort files
        if fileID(fastqFile)[1] is '1':
            if fileID(fastqFile)[2] == 'paired':
                paired1=path2fastq
                gate1=fileID(fastqFile)[0]
            elif fileID(fastqFile)[2] == 'unpaired':
                paired2=path2fastq
        if fileID(fastqFile)[1] is '2':
            if fileID(fastqFile)[2] == 'paired':
                paired3=path2fastq
            elif fileID(fastqFile)[2] == 'unpaired':
                paired4=path2fastq
                gate2=fileID(fastqFile)[0]
                if gate1 == gate2:
                    command='kallisto quant -i %s -o %s%s %s %s'%(genomeIndex,kallistoDir,gate1,paired1,paired3)
                    os.system(command)


        
        
