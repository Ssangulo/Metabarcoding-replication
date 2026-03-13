# =============================================================================
# 1_demultiplexing.R
# Demultiplexing and primer trimming of ITS amplicon sequencing data
#
# Tools required (install via conda):
#   - sabre  (barcode-based demultiplexing)
#   - cutadapt (primer/adapter trimming)
#   - dos2unix (line-ending conversion)
# =============================================================================

# ---- CONDA ENVIRONMENT SETUP ------------------------------------------------
# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:

# cutadapt
conda install -c bioconda cutadapt
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

# ---- RAW DATA SETUP ---------------------------------------------------------

#create data folder in my data directory
mkdir data
cd data 
mkdir seqfiles

# Open the terminal and navigate to your sequencing raw_data folder. For me it is here 
cd /data/bigexpansion/michadm/seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_/X204SC24050751-Z01-F001/01.RawData/
cd /data/bigexpansion/michadm/seqdata/2024-06-20_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding__batch2/X204SC24050751-Z01-F002/01.RawData/


# Copy raw data folder to my data directory 
cp -r /data/bigexpansion//seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_/X204SC24050751-Z01-F001/01.RawData /data/lastexpansion//data/

#Create folders for each plate and split samples by plate

cd /data/lastexpansion/_ang/data/rawdata/
mkdir P1 P2 P3 P4
cp -r *_P1 /data/lastexpansion/_ang/data/newrawdata01/P1/
cp -r *_P2 /data/lastexpansion/_ang/data/newrawdata01/P2/
cp -r *_P3 /data/lastexpansion/_ang/data/newrawdata01/P3/
cp -r *_P4 /data/lastexpansion/_ang/data/newrawdata01/P4/

# ---- FILE CONCATENATION (two delivery batches) ------------------------------
# Since I received the raw data in two separated batches, need to concatenate the files
# Concatenating the two forward files together and the two reverse files together
# This script is for P1 - adapt for other plates

cd /data/lastexpansion/_ang/data/newrawdata01/P1/
mkdir seqfiles
cp **/*.fq.gz /data/lastexpansion/_ang/data/newrawdata01/P1/seqfiles

cd /data/lastexpansion/_ang/data/newrawdata01/P1/seqfiles

# Get a list of unique sample prefixes (up to the FKDL part of naming)
# THIS WILL WORK ONLY FOR P1 -> CHANGE SCRIPT FOR OTHER PLATES

for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P1)_.*_L1_[12]\.fq\.gz/\1/' | sort | uniq); do

# Concatenate forward reads for each sample
cat ${sample}*_L1_1.fq.gz > ${sample}_concat_1.fq.gz

# Concatenate reverse reads for each sample
cat ${sample}*_L1_2.fq.gz > ${sample}_concat_2.fq.gz

#Verify the concatenation
echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"
done

# Concatenation for P2:
for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P2)_.*_[12]\.fq\.gz/\1/' | sort | uniq); do

cat ${sample}*_1.fq.gz > ${sample}_concat_1.fq.gz
cat ${sample}*_2.fq.gz > ${sample}_concat_2.fq.gz

echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"

done

# Move concatenated files to their destination
mv *_concat_1.fq.gz /data/lastexpansion/_ang/data/newrawdata01/P1/concatP1
mv *_concat_2.fq.gz /data/lastexpansion/_ang/data/newrawdata01/P1/concatP1

# ---- DEMULTIPLEXING WITH SABRE ----------------------------------------------
# Based on the forward barcode read for plate1 designated in the barcode file
# (euka01_repbarcodes.txt)

# If you made the barcode file with a text editor - use dos2unix to fix line endings
source ~/miniforge3/bin/activate 
conda activate myenv 

dos2unix euka01_repbarcodes.txt

echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /data/lastexpansion/_ang/data/euka01_repbarcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f.fq ${bn}_rep1f.fq;
mv rep1r.fq ${bn}_rep1r.fq;
mv rep2f.fq ${bn}_rep2f.fq;
mv rep2r.fq ${bn}_rep2r.fq;
mv rep3f.fq ${bn}_rep3f.fq;
mv rep3r.fq ${bn}_rep3r.fq;
mv rep4f.fq ${bn}_rep4f.fq;
mv rep4r.fq ${bn}_rep4r.fq;  done' > eukasabre.sh

# Run shell script
bash eukasabre.sh

# ---- PRIMER TRIMMING WITH CUTADAPT -----------------------------------------
# Doing replicate by replicate; reverse primer made an exact match to the
# corresponding reverse barcode (not removed by sabre).
# rc = reverse complement
#
# Base primers:
#   ITS4  Fwd: TCCTCCGCTTATTGATATGC
#   ITS86F Rv: GTGAATCATCGAATCTTTGAA

# --- Plate 1, Replicate 1 (reverse barcode: AGGAA) ---
# ITS4    (fwd):   TCCTCCGCTTATTGATATGC
# ITS86F  (rcRv):  TTCAAAGATTCGATGATTCACTTCCT
# ITS4    (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F  (Rv):    AGGAAGTGAATCATCGAATCTTTGAA

echo 'for i in *rep1f.fq; do bn=${i/rep1f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACTTCCT -A AGGAAGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep1out1.fq.gz --untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz ${bn}rep1f.fq ${bn}rep1r.fq; done' > ITSrep1.sh
bash ITSrep1.sh

# --- Plate 1, Replicate 2 (reverse barcode: GAGTGG) ---
# ITS4    (fwd):   TCCTCCGCTTATTGATATGC
# ITS86F  (rcRv):  TTCAAAGATTCGATGATTCACCCACTC
# ITS4    (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F  (rv):    GAGTGGGTGAATCATCGAATCTTTGAA

echo 'for i in *rep2f.fq; do bn=${i/rep2f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACCCACTC -A GAGTGGGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep2out1.fq.gz --untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz ${bn}rep2f.fq ${bn}rep2r.fq; done' > ITSrep2.sh
bash ITSrep2.sh

# --- Plate 1, Replicate 3 (reverse barcode: CCACGTC) ---
# ITS4    (fwd):   TCCTCCGCTTATTGATATGC
# ITS86F  (rcRv):  TTCAAAGATTCGATGATTCACGACGTGG
# ITS4    (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F  (Rv):    CCACGTCGTGAATCATCGAATCTTTGAA

echo 'for i in *rep3f.fq; do bn=${i/rep3f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGACGTGG -A CCACGTCGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep3out1.fq.gz --untrimmed-paired-output ${bn}.rep3out2.fq.gz -o ${bn}.rep3.trim1.fq.gz -p ${bn}.rep3.trim2.fq.gz ${bn}rep3f.fq ${bn}rep3r.fq; done' > ITSrep3.sh
bash ITSrep3.sh

# --- Plate 1, Replicate 4 (reverse barcode: TTCTCAGC) ---
# ITS4    (fwd):   TCCTCCGCTTATTGATATGC
# ITS86F  (rcRv):  TTCAAAGATTCGATGATTCACGCTGAGAA
# ITS4    (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F  (Rv):    TTCTCAGCGTGAATCATCGAATCTTTGAA

echo 'for i in *rep4f.fq; do bn=${i/rep4f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGCTGAGAA -A TTCTCAGCGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep4out1.fq.gz --untrimmed-paired-output ${bn}.rep4out2.fq.gz -o ${bn}.rep4.trim1.fq.gz -p ${bn}.rep4.trim2.fq.gz ${bn}rep4f.fq ${bn}rep4r.fq; done' > ITSrep4.sh
bash ITSrep4.sh

# ---- MOVE TRIMMED FILES TO DADA2 INPUT DIRECTORY ---------------------------
# Make a new folder called P1trimmed and move all demultiplexed, oligo-trimmed files there.
mkdir P1trimmed
mv *trim1.fq.gz /data/lastexpansion/_ang/data/trimmed/P1trimmed/
mv *trim2.fq.gz /data/lastexpansion/_ang/data/trimmed/P1trimmed/

# Repeat the above sabre + cutadapt steps for plates P2, P3, P4,
# adapting the barcode sequences accordingly.
# Then proceed to script: 2_DADA2_lulu.R
