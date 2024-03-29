# genomics_workshop
Materials for UNSA workshop, roughly following the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop)

#### Things to install for Windows users
- VirtualBox (https://www.virtualbox.org/wiki/Downloads)
- A Linux-based OS vdi (e.g., the Ubuntu vdi at https://sourceforge.net/projects/osboxes/files/v/vb/55-U-u/19.04/19.0464.7z/download).

The vdi may need to be decompressed with Unarchiver (https://theunarchiver.com/). You can follow directions for using the vdi with VirtualBox at: https://blogs.oracle.com/oswald/importing-a-vdi-in-virtualbox

Once open, use the password "osboxes.org" to log in. Then find the terminal window in the Virtual Desktop and open that.

Once in the virtual machine install:
- Phyluce: https://phyluce.readthedocs.io/en/latest/installation.html
(instructions at the link above, note that it requires conda)

The following instructions roughly follow the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop).

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
    gzip *.fastq

#### 8. Illumiprocessor requires a very specific naming scheme, so change the names as follows:
- SRR5678249_1.fastq.gz -> SRR5678249_A_L001_R1_001.fastq.gz
- SRR5678249_2.fastq.gz -> SRR5678249_A_L001_R2_002.fastq.gz
- SRR5678258_1.fastq.gz -> SRR5678258_A_L001_R1_001.fastq.gz
- SRR5678258_2.fastq.gz -> SRR5678258_A_L001_R2_002.fastq.gz
Most of the name don't matter (only the first and last parts that are separated by underscores "_").

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
SqCL matches contigs to our reference loci using blat. You may have to download blat if you don't already have it. Blat executables are available at: http://hgdownload.soe.ucsc.edu/admin/exe/. Then you must specify the location of the blat executable after the --blat flag as I did above.

#### 14. Output a pseudo-reference genome containing the best sequence between our two samples for each locus
    python make_PRG_mod.py --lineage l1 --file sample_map_for_SqCL.csv --dir ~/Desktop/genomics_workshop --adir ~/Desktop/genomics_workshop/3_velvet-output/contigs --mdir ~/Desktop/genomics_workshop/4_match-contigs-to-probes --keep easy_recip_match --outdir ~/Desktop/genomics_workshop/4_match-contigs-to-probes
If all we cared about was the consensus sequence for each individual, we could stop here! Since we want to call all variants and genotype each individual, we need to keep going!

#### 15. Now do some things to prepare the sequences for later mapping/variant-calling
Create a dictionary from the reference contigs
```
java -Xmx2g -jar ~/anaconda/pkgs/picard-1.106-0/jar/CreateSequenceDictionary.jar \
    R=./4_match-contigs-to-probes/l1.fasta \
    O=./4_match-contigs-to-probes/l1.dict
```

Index the reference 
	samtools faidx ./4_match-contigs-to-probes/l1.fasta

Assign the PRG as an index for bwa:
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
	
#### 19. Add read groups to keep track of individuals after combining pileups
	java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar I=./6_picard/Xiphorhynchus_obsoletus_AMNH12343-pairedreads_FM.bam O=./6_picard/Xiphorhynchus_obsoletus_AMNH12343_RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=TestXX RGLB=Lib1 RGID=Xiphorhynchus_obsoletus_AMNH12343 RGSM=Xiphorhynchus_obsoletus_AMNH12343 VALIDATION_STRINGENCY=LENIENT
	java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar I=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642-pairedreads_FM.bam O=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642_RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=TestXX RGLB=Lib1 RGID=Xiphorhynchus_obsoletus_LSUMNS35642 RGSM=Xiphorhynchus_obsoletus_LSUMNS35642 VALIDATION_STRINGENCY=LENIENT
	
#### 20. Mark PCR duplicates (these do not add info and complicate genotyping)
	java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar I=./6_picard/Xiphorhynchus_obsoletus_AMNH12343_RG.bam O=./6_picard/Xiphorhynchus_obsoletus_AMNH12343_MD.bam METRICS_FILE=./6_picard/Xiphorhynchus_obsoletus_AMNH12343.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false
	java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar I=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642_RG.bam O=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642_MD.bam METRICS_FILE=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false
	
#### 21. Merge the BAM pileups
```
java -Xmx2g -jar ~/anaconda/jar/MergeSamFiles.jar \
    SO=coordinate \
    AS=true \
    I=./6_picard/Xiphorhynchus_obsoletus_AMNH12343_MD.bam \
    I=./6_picard/Xiphorhynchus_obsoletus_LSUMNS35642_MD.bam \
    O=./7_merge-bams/All.bam 
```

#### 22. Index the merged bam file
	samtools index ./7_merge-bams/All.bam 

#### 23. Call indels
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./7_merge-bams/All.bam  \
    --minReadsAtLocus 7 \
    -o ./8_GATK/All.intervals
```

#### 24. Realign pileup around indels
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./7_merge-bams/All.bam  \
    -targetIntervals ./8_GATK/All.intervals \
    -LOD 3.0 \
    -o ./8_GATK/All_RI.bam
```

#### 25. Call SNPs
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./8_GATK/All_RI.bam \
    -gt_mode DISCOVERY \
    -o ./8_GATK/All_raw_SNPs.vcf \
    -ploidy 2 \
    -rf BadCigar
```
Note that you may want to output ALL (including invariant) sites for some analyses. Do this by adding "-output_mode EMIT_ALL_SITES".

#### 26. Annotate SNPs
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./8_GATK/All_RI.bam \
    -G StandardAnnotation \
    -V:variant,VCF ./8_GATK/All_raw_SNPs.vcf \
    -XA SnpEff \
    -o ./8_GATK/All_raw_SNPs_annotated.vcf \
    -rf BadCigar      
 ```
   
#### 27. Annotate indels
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./8_GATK/All_RI.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o ./8_GATK/Genus_species_SNPs_indels.vcf \
    -rf BadCigar         
```

#### 28. Mask indels
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -V ./8_GATK/All_raw_SNPs.vcf \
    --mask ./8_GATK/Genus_species_SNPs_indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "BadValidation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 5.0" \
    --filterName "LowVQCBD" \
    -o ./8_GATK/All_SNPs_no_indels.vcf  \
    -rf BadCigar
```

#### 29. Filter SNPs
	cat ./8_GATK/All_SNPs_no_indels.vcf | grep 'PASS\|^#' > ./8_GATK/All_SNPs_pass-only.vcf 

#### 30. Phase SNPs
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -I ./8_GATK/All_RI.bam \
    --variant ./8_GATK/All_SNPs_pass-only.vcf \
    -L ./8_GATK/All_SNPs_pass-only.vcf \
    -o ./8_GATK/All_SNPs_phased.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar
```

WE NOW HAVE A FINAL PHASED, FILTERED VCF FILE OF SNPS ("All_SNPs_phased.vcf")! This serves as the basis (with minor formatting tweaks) for the input into many programs like Structure, adegenet, EEMS, etc (see some examples in the "Converters" folder). This is where most folks will stop. However, to analyze whole sequences from each locus (for e.g. phylogenetics/phylogeography) we may want to produce multi-sequence alignments. That process is described below.

There are various options for producing sequences for multi-sequence alignment. One way is to use the contigs from Velvet, the pseudo-reference genome. Each contig sequence represents a haplotype for a locus (even though we know many species are diploid or polyploid). A good pipeline for doing this is available at the Phyluce site: https://phyluce.readthedocs.io/en/latest/tutorial-one.html. Another options is to use the GATK variant calls to assemble diploid sequences ("diplotypes") for each individual so we can make multi-species alignments with both sequences at each locus. I plan to implement a strategy for doing this using just the GATK output, which will use the function "EMIT_ALL_SITES" from GATK to obtain both variant and invariant bases in order to construct sequences. In the meantime, we can use a shortcut method, inserting the variable sites back into the consensus sequences from Velvet. Although this provides us with two sequences per locus per individual, it erroneously represents the completeness of the data set (individuals missing data from a locus will have that data added in according to the consensus sequence from the pseudo-reference genome). This is not optimal, but it illustrates the multi-sequence alignment part of the pipeline.

#### 32. Make VCFs for each individual

Individual 1:
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -T SelectVariants \
    --variant ./8_GATK/All_SNPs_phased.vcf \
    -o ./8_GATK/Xiphorhynchus_obsoletus_AMNH12343_SNPs_phased.vcf \
    -sn Xiphorhynchus_obsoletus_AMNH12343 \
    -rf BadCigar
```
Individual 2:
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -T SelectVariants \
    --variant ./8_GATK/All_SNPs_phased.vcf \
    -o ./8_GATK/Xiphorhynchus_obsoletus_LSUMNS35642_SNPs_phased.vcf \
    -sn Xiphorhynchus_obsoletus_LSUMNS35642 \
    -rf BadCigar
```
    
#### 32. For each individual, make SNP tables (the current input format for SNP data for multi-sequence alignment steps)

Individual 1:
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -V ./8_GATK/Xiphorhynchus_obsoletus_AMNH12343_SNPs_phased.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o ./9_SNP-tables/Xiphorhynchus_obsoletus_AMNH12343_SNPs_phased-table.txt \
    -rf BadCigar
```

Individual 2:
```
java -Xmx2g -jar ~/anaconda/pkgs/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R ./4_match-contigs-to-probes/l1.fasta \
    -V ./8_GATK/Xiphorhynchus_obsoletus_LSUMNS35642_SNPs_phased.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o ./9_SNP-tables/Xiphorhynchus_obsoletus_LSUMNS35642_SNPs_phased-table.txt \
    -rf BadCigar
```

#### 33. For each individual, Add phased SNPs to reference and optionally filter (seqcap_pop script)

Individual 1:
```
python add_phased_snps_to_seqs_filter.py \
	./4_match-contigs-to-probes/l1.fasta \
	./9_SNP-tables/Xiphorhynchus_obsoletus_AMNH12343_SNPs_phased-table.txt \
	./10_sequences/Xiphorhynchus_obsoletus_AMNH12343_sequences.txt \
	1
```
Individual 2:
```
python add_phased_snps_to_seqs_filter.py \
	./4_match-contigs-to-probes/l1.fasta \
	./9_SNP-tables/Xiphorhynchus_obsoletus_LSUMNS35642_SNPs_phased-table.txt \
	./10_sequences/Xiphorhynchus_obsoletus_LSUMNS35642_sequences.txt \
	1
```

#### 34. Collate sequences from all individuals into files by UCE (seqcap_pop script)
```
python collate_sample_fastas_GATK.py \
	./10_sequences/ \
	./11_fasta-parts/ \
	sequences.txt
```

#### 35. Do a final alignment to make sure these sequences line up (MAFFT)
```
python run_mafft.py \
	./11_fasta-parts/ \
	./12_raw-alignments/
```

#### 36. Convert MAFFT output to Phylip alignments
```
python process_mafft_alignments_GATK.py \
	./12_raw-alignments/ \
	./13_processed-phylip/
```

DONE! These alignments can now be used as input for phyogenetic programs to build gene trees. Those gene trees can be used as input for e.g. Astral, *BEAST, and other gene tree-species tree methods.

