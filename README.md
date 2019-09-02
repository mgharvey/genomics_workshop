# genomics_workshop
Materials for UNSA workshop

Instructions roughly follow the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop)

### 1.	Make a directory for this workshop (e.g. "genomics_workshop"), along with a bunch of subfolders within that directory
- 1_raw-reads
- 2_clean-reads
- 3_velvet-output
- 4_match-contigs-to-probes
- 5_mapping
- 6_picard
- 7_merge-bams
- 8_GATK
- 9_SNP-tables
- 10_sequences
- 11_fasta-parts
- 12_raw-alignments
- 13_processed-phylip

Here are commands that will do this:
    cd ~/Desktop/
    mkdir genomics_workshop 
    cd genomics_workshop
    mkdir 1_raw-reads 2_clean-reads 3_velvet-output 4_match-contigs-to-probes 5_mapping 6_picard 7_merge-bams 8_GATK 9_SNP-tables 10_sequences 11_fasta-parts 12_raw-alignments 13_processed-phylip 

### 2. Install PHYLUCE
https://phyluce.readthedocs.io/en/latest/installation.html

### 3.	Download data files for two samples
- https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5678249 (LSUMNS_35642)
- https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5678258 (AMNH_12343)

### 4.	Download the NCBI SRA toolkit 
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

### 5.	Convert data files to FASTQ format (one file for each mate pair)
    /Users/michaelharvey/Applications/sratoolkit.2.9.6-1-mac64/bin/fastq-dump --split-files ~/Downloads/SRR5678249
    /Users/michaelharvey/Applications/sratoolkit.2.9.6-1-mac64/bin/fastq-dump --split-files ~/Downloads/SRR5678258

### 6. Move FASTQ files to the working directory
    cd ~/Downloads
    cp ~/Downloads/SRR*.fastq ~/Desktop/genomics_workshop/1_raw-reads/

### 7. Illumiprocessor requires compressed (gzipped) files, so gzip them
    cd  ~/Desktop/genomics_workshop/1_raw-reads/
    gzip *

### 8. Illumiprocessor requires a very specific naming scheme, so change the names as follows:
- SRR5678249_1.fastq.gz -> SRR5678249_A_L001_R1_001.fastq.gz
- SRR5678249_2.fastq.gz -> SRR5678249_A_L001_R2_002.fastq.gz
- SRR5678258_1.fastq.gz -> SRR5678258_A_L001_R1_001.fastq.gz
- SRR5678258_2.fastq.gz -> SRR5678258_A_L001_R2_002.fastq.gz

### 8. Make an Illumiprocessor config file (or see example file "illumiprocessor.conf")

### 9. Run Illumiprocessor to clean the raw reads
    illumiprocessor --input ./1_raw-reads --output ./2_clean-reads/ --config illumiprocessor.conf 

### 10. Make a Velvet config file (simply containing a header, the sample name, and location of cleaned reads)

> [samples]
>
> Xiphorhynchus_obsoletus_AMNH12343:/Users/michaelharvey/Desktop/genomics_workshop/2_clean-reads/Xiphorhynchus_obsoletus_AMNH12343/split-adapter-quality-trimmed

### 11. Run Velvet to de novo assemble our reads into contigs
    phyluce_assembly_assemblo_velvet --conf velvet.conf --output 3_velvet-output --clean --cores 4 --kmer 67

WE NOW HAVE A REFERENCE ASSEMBLY! CHECK IT OUT!

### 12. Match our de novo contigs to the loci we targeted in our laboratory hybridization and enrichment
    phyluce_assembly_match_contigs_to_probes 

