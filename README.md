# genomics_workshop
Materials for UNSA workshop, roughly following the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop)

#### Things to install for Windows users
- VirtualBox (https://www.virtualbox.org/wiki/Downloads)
- An Ubuntu vdi (e.g., https://sourceforge.net/projects/osboxes/files/v/vb/55-U-u/19.04/19.0464.7z/download).

The vdi may need to be decompressed with Unarchiver (https://theunarchiver.com/). You can follow directions for using the vdi with VirtualBox at: https://blogs.oracle.com/oswald/importing-a-vdi-in-virtualbox

Once open, use the password "osboxes.org" to log in. Then find the terminal window in the Virtual Desktop and open that.

Once in the virtual machine install:
- Phyluce: https://phyluce.readthedocs.io/en/latest/installation.html
(requires conda)

Instructions roughly follow the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop)

#### 1.	Make a directory for this workshop (e.g. "genomics_workshop"), along with a bunch of subfolders within that directory
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
```
cd ~/Desktop/
mkdir genomics_workshop 
cd genomics_workshop
mkdir 1_raw-reads 2_clean-reads 3_velvet-output 4_match-contigs-to-probes 5_mapping 6_picard 7_merge-bams 8_GATK 9_SNP-tables 10_sequences 11_fasta-parts 12_raw-alignments 13_processed-phylip 
```

#### 2. Install PHYLUCE
https://phyluce.readthedocs.io/en/latest/installation.html

#### 3.	Download data files for two samples
- https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5678249 (LSUMNS_35642)
- https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5678258 (AMNH_12343)

#### 4.	Download the NCBI SRA toolkit 
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

#### 5.	Convert data files to FASTQ format (one file for each mate pair)
    /Users/michaelharvey/Applications/sratoolkit.2.9.6-1-mac64/bin/fastq-dump --split-files ~/Downloads/SRR5678249
    /Users/michaelharvey/Applications/sratoolkit.2.9.6-1-mac64/bin/fastq-dump --split-files ~/Downloads/SRR5678258

#### 6. Move FASTQ files to the working directory
    cd ~/Downloads
    cp ~/Downloads/SRR*.fastq ~/Desktop/genomics_workshop/1_raw-reads/

#### 7. Illumiprocessor requires compressed (gzipped) files, so gzip them
    cd  ~/Desktop/genomics_workshop/1_raw-reads/
    gzip *

#### 8. Illumiprocessor requires a very specific naming scheme, so change the names as follows:
- SRR5678249_1.fastq.gz -> SRR5678249_A_L001_R1_001.fastq.gz
- SRR5678249_2.fastq.gz -> SRR5678249_A_L001_R2_002.fastq.gz
- SRR5678258_1.fastq.gz -> SRR5678258_A_L001_R1_001.fastq.gz
- SRR5678258_2.fastq.gz -> SRR5678258_A_L001_R2_002.fastq.gz

#### 9. Make an Illumiprocessor config file (or see example file "illumiprocessor.conf")

#### 10. Run Illumiprocessor to clean the raw reads
    illumiprocessor --input ./1_raw-reads --output ./2_clean-reads/ --config illumiprocessor.conf 

#### 11. Make a Velvet config file (simply containing a header, the sample name, and location of cleaned reads)

> [samples]  
> Xiphorhynchus_obsoletus_AMNH12343:/Users/michaelharvey/Desktop/genomics_workshop/2_clean-reads/Xiphorhynchus_obsoletus_AMNH12343/split-adapter-quality-trimmed
> Xiphorhynchus_obsoletus_LSUMNS35642:/Users/michaelharvey/Desktop/genomics_workshop/2_clean-reads/Xiphorhynchus_obsoletus_LSUMNS35642/split-adapter-quality-trimmed

#### 12. Run Velvet to "de novo" assemble our reads into contigs for both samples
    phyluce_assembly_assemblo_velvet --conf velvet.conf --output 3_velvet-output --clean --cores 4 --kmer 67

WE NOW HAVE DE NOVO ASSEMBLIES! Open the contigs.fa files in 3_velvet-output and check them out.

#### 13. Match our de novo contigs to the loci we targeted in our laboratory hybridization and enrichment - the scripts we use for this step and the next one are in the Github directory and are modified from Sonal Singhal's excellent SqCL pipeline (https://github.com/singhal/SqCL).
    python match_contigs_to_probes_mod.py --blat /usr/local/bin/blat --sample Xiphorhynchus_obsoletus_AMNH12343 --dir ~/Desktop/genomics_workshop --evalue 1e-30 --db uce-and-exon-probes.fasta --outdir ./4_match-contigs-to-probes/
    python match_contigs_to_probes_mod.py --blat /usr/local/bin/blat --sample Xiphorhynchus_obsoletus_LSUMNS35642 --dir ~/Desktop/genomics_workshop --evalue 1e-30 --db uce-and-exon-probes.fasta --outdir ./4_match-contigs-to-probes/
SqCL matches contigs to our reference loci using blat. We may have to install blat if it's not installed already. Blat executables are available at: http://hgdownload.soe.ucsc.edu/admin/exe/

#### 14. Output a pseudo-reference genome containing the best sequence between our two samples for each locus
    python make_PRG_mod.py --lineage l1 --file sample_map_for_SqCL.csv --dir ~/Desktop/genomics_workshop --adir ~/Desktop/genomics_workshop/3_velvet-output/contigs --mdir ~/Desktop/genomics_workshop/4_match-contigs-to-probes --keep easy_recip_match --outdir ~/Desktop/genomics_workshop/4_match-contigs-to-probes

#### 15. Assign a reference sequence for bwa to use
    bwa index -a is ./4_match-contigs-to-probes/l1.fasta
    
#### 16. Map reads from each individual to the reference using bwa (this makes a pileup)

Individual 1:
```
bwa mem -t 4 ./4_match-contigs-to-probes/l1.fasta /Users/michaelharvey/Desktop/genomics_workshop/2_clean-reads/Xiphorhynchus_obsoletus_AMNH12343/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_AMNH12343-READ1.fastq.gz /Users/michaelharvey/Desktop/genomics_workshop/2_clean-reads/Xiphorhynchus_obsoletus_AMNH12343/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_AMNH12343-READ2.fastq.gz > /Users/michaelharvey/Desktop/genomics_workshop/5_mapping/Xiphorhynchus_obsoletus_AMNH12343-pairedreads.sam
bwa mem -t 4 ./4_match-contigs-to-probes/l1.fasta ./2_clean-reads/Xiphorhynchus_obsoletus_AMNH12343/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_AMNH12343-READ-singleton.fastq.gz > ./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-unpairedreads.sam
```

Individual 2:
```
bwa mem -t 4 ./4_match-contigs-to-probes/l1.fasta ./2_clean-reads/Xiphorhynchus_obsoletus_LSUMNS35642/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_LSUMNS35642-READ1.fastq.gz ./2_clean-reads/Xiphorhynchus_obsoletus_LSUMNS35642/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_LSUMNS35642-READ2.fastq.gz > ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads.sam
bwa mem -t 4 ./4_match-contigs-to-probes/l1.fasta ./2_clean-reads/Xiphorhynchus_obsoletus_LSUMNS35642/split-adapter-quality-trimmed/Xiphorhynchus_obsoletus_LSUMNS35642-READ-singleton.fastq.gz > ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-unpairedreads.sam
```

#### 17. Convert the SAM pileup to BAM format
```
samtools view -bS ./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-pairedreads.sam > ./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-pairedreads.bam
samtools view -bS ./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-unpairedreads.sam > ./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-unpairedreads.bam

samtools view -bS ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads.sam > ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads.bam
samtools view -bS ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-unpairedreads.sam > ./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-unpairedreads.bam
```

#### 18. Fix mate pair information
	java -Xmx2g -jar ~/anaconda/jar/FixMateInformation.jar I=./5_mapping/Xiphorhynchus_obsoletus_AMNH12343-pairedreads.bam O=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-pairedreads_FM.bam 
	java -Xmx2g -jar ~/anaconda/jar/FixMateInformation.jar I=./5_mapping/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads.bam O=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads_FM.bam 
	
#### 19. Sort
	samtools sort -@ 4 -o ./6_picard/Xiphorhynchus_obsoletus_AMNH12343-pairedreads_ST.bam ./6_picard/Xiphorhynchus_obsoletus_AMNH12343-pairedreads_FM.bam
	samtools sort -@ 4 -o ./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads_ST.bam ./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads_FM.bam

#### 20. Add read groups to keep track of individuals after combining pileups
	java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar I=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-aln_FM.bam O=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-aln_RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=TestXX RGLB=Lib1 RGID=Xiphorhynchus_obsoletus_AMNH12343 RGSM=Xiphorhynchus_obsoletus_AMNH12343 VALIDATION_STRINGENCY=LENIENT
	java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar I=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-aln_FM.bam O=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-aln_RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=TestXX RGLB=Lib1 RGID=Xiphorhynchus_obsoletus_LSUMNS35642 RGSM=Xiphorhynchus_obsoletus_LSUMNS35642 VALIDATION_STRINGENCY=LENIENT
	
#### 21. Mark PCR duplicates (these do not add info and complicate genotyping)
	java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar I=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-aln_RG.bam O=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-aln_MD.bam METRICS_FILE=./6_picard/Xiphorhynchus_obsoletus_AMNH12343.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false
	java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar I=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-aln_RG.bam O=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-aln_MD.bam METRICS_FILE=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false
	
#### 22. Create a dictionary from the reference contigs
	java -Xmx2g -jar ~/anaconda/pkgs/picard-1.106-0/jar/CreateSequenceDictionary.jar \
    R=./4_match-contigs-to-probes/PRG/l1.fasta \
    O=./4_match-contigs-to-probes/PRG/l1.dict

#### 23. Index the reference
	samtools faidx ./4_match-contigs-to-probes/PRG/l1.fasta

#### 24. Index the current pileup
	samtools index ./6_picard/Xiphorhynchus_obsoletus_AMNH12343-aln_MD.bam
	samtools index ./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-aln_MD.bam
