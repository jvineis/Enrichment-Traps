library('dada2');packageVersion('dada2')

##REMEMBER TO CHANGE THIS DEPENDING ON THE ANALYSIS YOU ARE CONDUCTING:
#path = "~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/"
path = "/Users/joevineis/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS" # change to this in order to process the other run - it has a very different error profile and many reads are lost if this value is not changed.


#sort files into R1 (forward) and R2 (reverse)
fnFs <- sort(list.files(path, pattern="-R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="-R2.fastq", full.names = TRUE))

#sample names
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

#quality plots
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


**** OUTPUT OF THE ABOVE PIPELINE
#input filtered denoisedF denoisedR merged nonchim
#J3-3-control 36313    35145     34345     34575  27631    9487
#J3-3-Red     38858    37508     35930     36803  27930   10851

#Create a table to save as an rds file that will be used to merge the two runs
seqtab = makeSequenceTable(mergers)

##REMEMBER TO CHANGE THIS DEPENDING ON THE ANALYSIS YOU ARE CONDUCTING:
saveRDS(seqtab, "~/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS/seqtab-comp-dilution.rds")
#saveRDS(seqtab, "~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/seqtab-original.rds")


##REMEMBER TO CHANGE THIS DEPENDING ON THE ANALYSIS YOU ARE CONDUCTING:

#write.table(seqtab.nochim, "~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/x_original-cultivation-seqtab.nochim.txt", sep = "\t")
write.table(seqtab.nochim, "~/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS/x_comp-dilution-seqtab.nochim.txt", sep = "\t")


**** HERE IS HOW TO MERGE THE TWO RUNS
st1 = readRDS("~/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS/seqtab-comp-dilution.rds")
st2 = readRDS("~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/seqtab-original.rds")

st.all = mergeSequenceTables(st1, st2)

seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)


write.table(seqtab, "~/Dropbox/TRAP-COMPETITION/seqtab-original-comp-dil-combined-runs.txt", sep = "\t")

#2. Then you need to fix the output using some type of txt editor so that the
# rows are sequences and the colums are samples..   Its just a matter of transpoosing the matrix and adding the "sample"
# as the column1 header. like this.

sample 15cm_S_Alt_2_1  15cm_S_Alt_2_2  15cm_S_Alt_2_3
ACTTGA 0  26  33
ACTTAA 0  100 30
AACTGG 25 33  26

# Save it as a text file called dada2-t-count-table.txt
## then use this handy script
# python ~/scripts/create-fasta-from-seqs.py -n dada2-t-count-table.txt -fa dada2-fasta.fa -o DADA2-MATRIX.txt
# which produces a count table with simple names instead of the asv as the identifier and a fasta file that you can use to identify the taxonomy
## I ran vsearch on this file to generate taxonmy like this in the terminal
#vsearch --usearch_global dada2-fasta.fa --db ~/scripts/databas/silva119.fa --blast6out NODE-HITS.txt --id 0.6
## Then I ran the wonderful script below in the terminal
#python ~/scripts/mu-dada2-phyloseq-creator.py -hits NODE-HITS.txt -tax_ref ~/scripts/databas/silva_fix.tax -dada2 DADA2-MATRIX.txt -fa dada2-fasta.fa
## which creates the PHYLOSEQ-TAX.txt and PHYLOSEQ-MATRIX.txt that you will use below
library(phyloseq)

mat = read.table("~/Dropbox/TRAP-COMPETITION/PHYLOSEQ-MATRIX.txt", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("~/Dropbox/TRAP-COMPETITION/PHYLOSEQ-TAX.txt", header = TRUE, sep = ";", row.names = 1)

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)

physeq = phyloseq(OTU,TAX)

per = transform_sample_counts(physeq, function (x) x/sum(x)*100)

lowabundnames = filter_taxa(per, function(x) mean(x) > 0.1) ## This filters out anything that has less than a mean of 0.2% relative abundance acrocss all samples
per_abund = prune_taxa(lowabundnames, per)

## DO any other filtering that you think should be required here.  

# Write the taxa to a file
write.table(tax_table(per_abund), "~/Dropbox/TRAP-COMPETITION/DADA2-TAX-meangreaterthanpoint1.txt", sep = '\t')
# Write the matrix
write.table(otu_table(per_abund), "~/Dropbox/TRAP-COMPETITION/DADA2-MATRIX-phyloseq-meangreaterthanpoint1.txt", sep = '\t')

### READ in a table that contains the percent relative abundance of the crosses experiment and export a tree
library(ape)
library(vegan)

dat = read.table("~/Dropbox/TRAP-COMPETITION/CROSSES-perabund-table.txt", header = TRUE, row.names = 1)
dist = vegdist(dat, method = "bray")
dclust = hclust(dist)
write.tree(as.phylo(dclust), "~/Dropbox/TRAP-COMPETITION/CROSSES-perabund-table.tre")

## READ in a table that contains the percente relative abundance of ALL samples and export a tree
dat = read.table("~/Dropbox/TRAP-COMPETITION/ALL-SAMPLES-perabund-table.txt", header = TRUE, row.names = 1)
dist = vegdist(dat, method = "bray")
dclust = hclust(dist)
write.tree(as.phylo(dclust), "~/Dropbox/TRAP-COMPETITION/ALL-SAMPLES-perabund-table.tre")


