##########################
# Ensuite on utilise Vsearch pour construire la table des OTUs
# https://grunwaldlab.github.io/rps10_barcode/03--abundance_matrix_preparation.html

library(Biostrings)
library(dada2)
library(parallel)
#install.packages("metacoder")
library(metacoder)
library(purrr)
library(brio)
library(stringr)
library(dplyr)
set.seed(100) #

source("~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo/1_dada2_illu_Ivan_mai2025.R")

setwd("~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo")
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# plot(table(nchar(getSequences(seqtab))))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab) # % de reads conservés (non-chimériques)
#plot(table(nchar(getSequences(seqtab.nochim))))

write.csv(t(seqtab.nochim),"/Users/tonyrobinet/sync/BioLav/M2_Yvan/Stage_M2/bioinfo/seqtabnochim_Amph16S.csv") # change with your path
uniquesToFasta(seqtab.nochim,"/Users/tonyrobinet/sync/BioLav/M2_Yvan/Stage_M2/bioinfo/seqtabnochim_Amph16S.fasta") # change with your path

write.table(t(seqtab.nochim), "/Users/tonyrobinet/sync/BioLav/M2_Yvan/Stage_M2/bioinfo/seqtabnochim_Amph16S.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) # change with your path
uniquesToFasta(seqtab.nochim, "/Users/tonyrobinet/sync/BioLav/M2_Yvan/Stage_M2/bioinfo/seqtabnochim_Amph16S.fna", ids=colnames(seqtab.nochim)) # change with your path


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#write.table(track,"~/alice/seqtabnochim_mangrove_16S_18S_ITS_track.txt") # change with your path


seq_abund <- readDNAStringSet("./seqtabnochim_Amph16S.fasta")
seqs <- read.table("./seqtabnochim_Amph16S.fna")
NCBI_ENA_Mito_DB <- readDNAStringSet("./RefBase_Amph16S/ena_ncbi/NCBI_ENA_Amphibiens_BZH_ET_invasifs_CROAA__1070F_1340R_assignTaxonomy.fasta.gz")



merged_read_seqs <- unlist(map(mergers, function(x) {
  x$sequence[x$accept]
})) # avec l'obj 'mergers' issu de mergePairs() du script 1_dada2

unique_merged_read_seqs <- unique(merged_read_seqs)
length(unique_merged_read_seqs)

raw_abundance_data <- makeSequenceTable(mergers)
hist(nchar(getSequences(raw_abundance_data)))
rowSums(raw_abundance_data)

## Create OTU abundance matrix
its_clustering_threshold <- 0.97

vsearch_cluster <- function(seqs, seq_abund, id_threshold = its_clustering_threshold, method = "fast") {
  # Check that VSEARCH is installed
  tryCatch(system2("vsearch", args = "--version", stdout = FALSE, stderr = FALSE),
           warning=function(w) {
             stop("vsearch cannot be found on PATH. Is it installed?")
           })
  
  # Run VSEARCH
  # seqs <- seqs[order(seq_abund, decreasing = TRUE)]
  input_fasta_path <- tempfile()
  write_lines(paste0('>', seq_along(seqs), ';size=', seq_abund, '\n', seqs), path = input_fasta_path)
  otu_centroid_path <- tempfile()
  command_args <- paste(paste0("--cluster_", method), 
                        input_fasta_path,
                        "--threads", detectCores() - 1,
                        "--id", id_threshold,
                        "--sizein",
                        "--strand plus",
                        "--fasta_width 0", # 0 = no wrapping in fasta file
                        "--centroids", otu_centroid_path)
  system2("vsearch", args = command_args, stdout = FALSE, stderr = FALSE)
  
  # Return OTU sequences
  centroids <- read_fasta(otu_centroid_path)
  names(centroids) <- str_match(names(centroids), pattern = 'size=(.+)$')[, 2]
  return(centroids)
}


unique_read_counts <- map_dbl(unique_merged_read_seqs, function(s) {
  sum(map_dbl(mergers, function(sample_data) {
    sum(sample_data$abundance[sample_data$sequence == s & sample_data$accept])
  }))
})

otu_seqs <- vsearch_cluster(seqs = unique_merged_read_seqs,
                                seq_abund = unique_read_counts,
                                id_threshold = its_clustering_threshold,
                                method = 'size') %>% toupper()


otus_per_sample <- map(rownames(raw_abundance_data), function(sample) {
  sample_id <- str_match(sample, pattern = '^(.+)_.+$')[, 2]
  merged_read_data <- mergers[[sample]]
  sample_otu_counts <- map_int(otu_seqs, function(s) {
    sum(merged_read_data$abundance[merged_read_data$sequence == s & merged_read_data$accept])
  })
  names(sample_otu_counts) <- otu_seqs
  all_unique_otus <- unique(otu_seqs)
  out <- as.integer(rep(0, length(all_unique_otus)))
  names(out) <- all_unique_otus
  out[names(sample_otu_counts)] <- sample_otu_counts
  out
  return(out)
})



raw_otu_abundance_data <- do.call(rbind, otus_per_sample)
rownames(raw_otu_abundance_data) <- rownames(raw_abundance_data) %>% substr(21,22)
rowSums(raw_otu_abundance_data)
length(raw_otu_abundance_data)


####################################################
## Assign taxonomy
# après avoir mis en forme le fichier de séquences références, cf. "creation_base_NCBI.ipynb"

setwd("~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo")

# Assignation to genus
tax_results_otu_genus <- assignTaxonomy(raw_otu_abundance_data, 
                                      refFasta = "./RefBase_Amph16S/ena_ncbi/NCBI_ENA_Amphibiens_BZH_ET_invasifs_CROAA__1070F_1340R_assignTaxonomy.fasta.gz", 
                                      taxLevels = c("Domaine", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                      minBoot = 50,
                                      tryRC = TRUE,
                                      outputBootstraps = TRUE,
                                      multithread = TRUE)


tax_results_otu_genus_plus <- addSpecies(tax_results_otu_genus[[1]], 
                                         "~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo/RefBase_Amph16S/ena_ncbi/NCBI_ENA_Amphibiens_BZH_ET_invasifs_CROAA__1070F_1340R_assignSpecies.fasta.gz",
                                         verbose=TRUE)
unname(tax_results_otu_genus_plus)


assignTax_as_char <- function(res) {
  out <- vapply(1:nrow(res$tax), FUN.VALUE = character(1), function(i) {
    paste(res$tax[i, ],
          res$boot[i, ],
          colnames(res$tax), 
          sep = '--', collapse = ';')
  })
  names(out) <- rownames(res$tax)
  return(out)
}


seq_tax_otu_genus <- c(assignTax_as_char(tax_results_otu_genus), assignTax_as_char(tax_results_otu_genus))

formatted_abund_otu <- t(raw_otu_abundance_data)
colnames(formatted_abund_otu) <- sub(colnames(formatted_abund_otu), pattern = "_.+$", replacement = "")

#formatted_abund_otu <- cbind(sequence = rownames(formatted_abund_otu), 
#                             taxonomy = tax_results_otu_genus_plus[rownames(formatted_abund_otu)], 
#                             formatted_abund_otu)

formatted_abund_otu <- cbind(sequence = rownames(formatted_abund_otu), 
                             taxonomy = cbind(formatted_abund_otu[,1],tax_results_otu_genus_plus), 
                             formatted_abund_otu)

formatted_abund_otu <- as_tibble(formatted_abund_otu)
formatted_abund_otu[2:10,2:10]

write.csv(formatted_abund_otu, "./formatted_abund_otu_97.csv")


###### Assignation to species

seq_uniques <- getUniques(seqtab.nochim)
path2 <- "~/sync/BioLav/M2_Ivan/Stage_M2/bioinfo/RefBase_Amph16S/ena_ncbi/NCBI_ENA_Amphibiens_BZH_ET_invasifs_CROAA__1070F_1340R_assignSpecies.fasta.gz"

tax_results_otu_species <- assignSpecies(seq_uniques, 
                                         refFasta = path2,
                                         allowMultiple = FALSE,
                                         tryRC = FALSE,
                                         n = 100,
                                         verbose = TRUE)


unname(tax_results_otu_species)








