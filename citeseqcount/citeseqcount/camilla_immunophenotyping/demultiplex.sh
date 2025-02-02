########## MBU_zjg20_covid ##########

### do demultiplexing

# see: https://genomicsequencing.cruk.cam.ac.uk/sequencing/help/Help!demultiplexing.action

# input file: R1 and R2 read of GEX only (this is where the index is)
# MTOX_07 TotalSeq A0451 anti-Nuclear Pore Complex Proteins Hashtag 1. Cat number 682205. Barcode Sequence: TTCCTGCCATTACTA

# first, create the index.csv file (as in link)


#### DETAILS

# experimental TotalSeqA kit: HTO B7



#### set up environment

DEMULTIPLEXER_TOOL=/suffolk/RawGenomicsH/Genomics_Group/software/demultiplexer.rhel/demuxFQ
WORKING_DIRECTORY=/suffolk/WorkGenomicsD/ec474/demultiplex/
INPUT_GEX_FASTQ=/suffolk/RawGenomicsH/Genomics_Group/fastqs/ImmunoPheno_project/007_Immuno_GP0030/SLX-22953_GEX/

CELL_RANGER_COUNT=/suffolk/RawGenomicsH/Genomics_Group/software/cellranger-7.1.0/cellranger
REFERENCES_MASKED=/suffolk/RawGenomicsH/Genomics_Group/references/masked/



########## STEP 1: DEMULTIPLEX ##########

# R1
${DEMULTIPLEXER_TOOL} -c -d -i -e -t 1 -r 0.01 -o ${WORKING_DIRECTORY} \
	-b SLX-22953.SITTC1.HLNNMDRX3.lostReads.R1.afterDemultiplexing.fq.gz \
	-s SLX-22953_demultiplexsummary_R1.txt \
	/suffolk/WorkGenomicsF/cl873/MBU_zjg20_covid/demultiplex/MBU_zjg20_covid_demultiplex_index_R1.tsv \
	${INPUT_GEX_FASTQ}SLX-22953.HLNNMDRX3.s_1.r_1.lostreads.fq.gz


# R2
${DEMULTIPLEXER_TOOL} -c -d -i -e -t 1 -r 0.01 -o ${WORKING_DIRECTORY} \
	-b SLX-22953.SITTC1.HLNNMDRX3.lostReads.R2.afterDemultiplexing.fq.gz \
	-s SLX-22953_demultiplexsummary_R2.txt \
	${WORKING_DIRECTORY}MBU_zjg20_covid_demultiplex_index_R2.tsv \
	${INPUT_GEX_FASTQ}SLX-22953.HLNNMDRX3.s_1.r_2.lostreads.fq.gz





########## STEP 2: CELL RANGER ##########

# see: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#ab-totalseqa

# first create 2 files: library.csv and feature_ref.csv

# then align to genome with CellRanger

# NUMT masked genome only



${CELL_RANGER_COUNT} count --id=SLX-22953_demultiplexed_masked \
			   --libraries=${WORKING_DIRECTORY}libraries.csv \
			   --transcriptome=${REFERENCES_MASKED}hg38/GRCh38 \
			   --feature-ref=${WORKING_DIRECTORY}feature_ref.csv \
			   --chemistry=ARC-v1




