# this a script to convert sra files to fastq files
import os,sys

# 0. user defined variables
sraDir='/Users/user/athan/transcriptomics/data/sra/'
fastqDir='/Users/user/athan/transcriptomics/data/fastq/'

# 1. select input files
sraFiles=os.listdir(sraDir)

# 2. define samples as single or paired end
singleEndSamples=['81','83','85','87']

# 3. execute fastq dump command
for sraFile in sraFiles:
    path2sra=sraDir+sraFile

    lastNumbers=sraFile.split('.')[0][-2:]
    if lastNumbers in singleEndSamples:
        command='fastq-dump --outdir %s --gzip --skip-technical  --readids --dumpbase --clip %s'%(fastqDir,path2sra)
    else:
        command='fastq-dump --outdir %s --gzip --skip-technical  --readids --dumpbase --split-files --clip %s'%(fastqDir,path2sra)

    os.system(command)
    print(sraFile + ' has been converted')
    sys.exit()
