en # Script to clean reads using trimmomatic
import os,sys

# User defined variables
fastqDir='/Users/user/athan/transcriptomics/data/fastq/'
trimJar='/Users/user/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar'
cleanDir='/Users/user/athan/transcriptomics/data/cleanfastq/'

# Obtain input files
fastqFiles=os.listdir(fastqDir)

# Identify single end samples
singleEndSamples=['81', '83', '85', '87']

# Run Trimmomatic
for fastqFile in fastqFiles:
    path2fastq=fastqDir+fastqFile
    path2clean=cleanDir+fastqFile
    fileName=fastqFile.split('.')[0]
    
    # Gather last numbers
    lastNumbers=fastqFile.split('.')[0][7:]
    # Separate single and paired end samples
    if lastNumbers in singleEndSamples:
        command='java -jar %s SE -phred33 %s %s ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(trimJar,path2fastq,path2clean)
        os.system(command)
    else:
         Identify forward and backward reads
        pairedNumber=lastNumbers.split('_')[1][0:]
        if pairedNumber is '1':
            forwardFile=fastqFile
            forwardFileName=fileName
        else:
            backwardFile=fastqFile
            backwardFileName=fileName

            path2fastq_1=fastqDir+forwardFile
            path2fastq_2=fastqDir+backwardFile
            path2clean_1=cleanDir+forwardFileName
            path2clean_2=cleanDir+backwardFileName
    
            command='java -jar %s PE -phred33 %s %s %s_paired.fastq.gz %s_unpaired.fastq.gz %s_paired.fastq.gz %s_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(trimJar,path2fastq_1,path2fastq_2,path2clean_1,path2clean_1,path2clean_2,path2clean_2)
            os.system(command)

