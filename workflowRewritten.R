
# CHANGE to the directory containing the fastq files after unzipping.
miseq_path <- "./MiSeq_SOP" 
list.files(miseq_path)


fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]
fnRs[1:3]

plotQualityProfile(fnFs[1:2]) #The first two forward reads
plotQualityProfile(fnRs[1:2]) #The first two reverse reads

#We define the filenames for the filtered fastq.gz files:
 
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


#Filter the forward and reverse reads**:

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


#Infer sequence variants
#we use the high-resolution DADA2 method to to
#infer amplicon sequence variants (ASVs) exactly, without imposing any
#arbitrary threshhold, and thereby resolving variants that differ by as
#little as one nucleotide

### Dereplication 

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# learning errors to distinguish between sequencing errors and biological variants

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF)
plotErrors(errR)

#In order to verify that the error rates have been reasonably
#well-estimated, we inspect the fit between the observed error rates
#(black points) and the fitted error rates (black lines) in. These figures show the frequencies of
#each type of transition as a function of the quality.

###########################SO CHECK THE PLOTS! ####################################

#The DADA2 sequence inference method can run in two different modes:
#Independent inference by sample (pool=FALSE), and inference from the
#pooled sequencing reads from all samples (pool=TRUE). Independent
#inference has the advantage that computation time is linear in the
#number of samples, and memory requirements are flat with the number of
#samples. This allows scaling out to datasets of almost unlimited size.

#Pooled inference is more computationally taxing, however, pooling
#improves the detection of --> rare variants <-- that were seen just once or
#twice in an individual sample but many times across all samples. As this
#dataset is not particularly large, we perform pooled inference.

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)

#Inspecting the dada-class object returned by dada:
#dadaFs[[1]]
i <- 1
for (i length(dadaFs)){
  print(length(dadaFs[[i]]$sequence), " real sequence variants")
  print(length(dadaFs[[i]]$map) "unique sequences in the",i , "th sample")
}

#The DADA2 sequence inference step removed (nearly) all substitution and
#indel errors from the data. We now merge together the inferred
#forward and reverse sequences, removing paired sequences that do not
#perfectly overlap as a final control against residual errors.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

#Construct sequence table and remove chimeras

seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
#table(nchar(getSequences(seqtabAll)))

#The error model did not address chimeric sequences, hence we expect this sequence table to include many chimeric sequences. 
#We now remove chimeric sequences by comparing each inferred sequence to the others in the table, and removing those that can be reproduced by stitching together two more abundant sequences.
seqtabNoC <- removeBimeraDenovo(seqtabAll)

#here chimera make up about
print('\nhere chimera make up about ', round(100*(ncol(seqtabAll)-ncol(seqtabNoC))/ncol(seqtabAll)),'% of the inferred sequence variants, but those variants account for only about',
      round(100*(sum(seqtabAll)-sum(seqtabNoC))/sum(seqtabAll)),'% of the total sequence reads')

#classification
#naive Bayesian classifier method for this purpose comparing sequence variants to a
#training set of classified sequences (here RDP v16)

#Additionally the dada2 tutorial website (https://benjjneb.github.io/dada2/training.html) contains formatted training fastas for:
#                                      - RDP training set, 
#                                     -  GreenGenes clustered at 97% identity, 
#                                      - Silva reference database] available.
#To follow this workflow, download the rdp_train_set_16.fa.gz file, and place it in the directory with the fastq files.

fastaRef <- "./rdp_train_set_16.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
#unname(head(taxTab))

#Construct phylogenetic tree
#We begin by performing a multiple-alignment using the Biocpkg "DECIPHER"
library(DECIPHER)
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#the pkg("phangorn")` is then used to construct a phylogenetic tree. 
#Here we first construct a neighbor-joining tree, and then fit a GTR+G+I 
#(Generalized time-reversible with Gamma rate variation)
#maximum likelihood tree using the neighbor-joining tree as a starting point.
library(phangorn)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



