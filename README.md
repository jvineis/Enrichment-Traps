# Enrichment-Traps: Analysis of 16S amplicon sequences and recovery of draft genomes from metagenomic data collected from microbial enrichment traps. 

## 16S Amplicon analysis




## Metagenomic analysis: assembly, mapping, anvio profiling, and draft genome recovery from each of the eight trap wells that we conducted metagenomic sequencing

### We begin with the assembly of quality filtered reads that can be obtained from NCBI under the submission id XXXX.  We assembled each of the metagenomic datasets recovered from each trap well. Below is an example of the assembly of the metagenomic data synthesized for trap well 1. We conducted the same assembly process for each of the eight metagenomic data sets  

    spades.py -k 21,33,55,77,99,127 --careful -1 cultivar-idx1-QUALITY_PASSED_R1.fastq -2 cultivar-idx1-QUALITY_PASSED_R2.fastq -o cultivar-idx1-spades

### Map each of the short read datasets to each individual assembly using bowtie2. Start by building a bowtie2 database and then run a for loop to map each of the short read datasets to the bowtie2 database.

    cd cultivar-idx1-spades
    
    bowtie2-build scaffolds.fasta bowtie2db
    
    for i in `cat samples.txt`; do clusterize bowtie2 --very-sensitive -x bowtie2db -1 ../"$i"-QUALITY_PASSED_R1.fastq -2 ../"$i"-QUALITY_PASSED_R2.fastq -S "$i".sam; done
    
### Build a contigs database for anvio using the scaffolds.fa file generated by the spades assembler

    anvi-gen-contigs-database -f scaffolds.fa -o contigs.db
    anvi-run-hmms -c contigs.db
    
### Profile each of the contigs.db using anvio. The samples.txt file is found within the data maintained here in the git.

    for i in `cat samples.txt`; do bowtie2 --very-sensitive -x bowtie2db -1 ../"$i"-QUALITY_PASSED_R1.fastq -2 ../"$i"-QUALITY_PASSED_R2.fastq -S "$i".sam; done
    
    for i in `cat samples.txt`; do samtools view -bS -F 4 "$i".sam > "$i".bam ; done
    
    for i in `cat samples.txt`; do anvi-init-bam "$i".bam -o "$i"-sorted.bam; done
    
    for i in `cat samples.txt`; do clusterize anvi-profile -c contigs.db -i "$i"-sorted.bam -o "$i"-PROFILE; done
    
### Merge the data from each assembly in order to create a profile db.  This profile and other required files to view the anvio display and the scaffold selections that we made for each draft genome. 

    anvi-merge *PROFILE/PROFILE.db -c contigs.db -o IDX1-MERGED 

### The files that you have generated using the above steps have produced all the data that you need to run the interactive display. These files are also available here https://figshare.com/account/home#/projects/97361 if you didn't make it through all the steps. However, the files found in figshare, contain the bin selections. Here is how to load the interactive display

    anvi-interactive -c contigs.db -p IDX1-MERGED/PROFILE.db
