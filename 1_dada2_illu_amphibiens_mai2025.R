# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("dada2")
library(dada2)

# On traite les trois  marqueurs (16SV4 18sV9 et ITS2) ensemble, on sort des seqtab_nochim communes
# puis on sépare les marqueurs dans QIIME2 avec $ qiime feature-classifier extract-reads
# Unknown_CP224-003M0001_good_2.fastq.gz

path <- "~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo/Amph16S_raw_filter_trim/" # change with your path

list.files(path)
setwd(path)
# path <- "/Volumes/BOUGREDANE/191224_mang1"
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_good_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_good_2.fastq.gz", full.names = TRUE))
fnFs <- sort(fnFs)
fnRs <- sort(fnRs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_good_1"), `[`, 1)
sample.names
#plotQualityProfile(fnFs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=100, matchIDs=TRUE,
                     maxN=0, maxEE=c(2,2), rm.phix=TRUE, #trimLeft=15, trimRight=15,
                     compress=TRUE, multithread=TRUE)
out

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
 dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

