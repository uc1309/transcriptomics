import os,sys

fastqDir='/Users/user/athan/transcriptomics/data/fastq/'

fastqFiles=os.listdir(fastqDir)

for fastqFile in fastqFiles:
    path2fastq=fastqDir+fastqFile

    command = 'gunzip %s'%(path2fastq)
    os.system(command)
    
