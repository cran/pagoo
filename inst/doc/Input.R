## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE,  warning=FALSE, message=FALSE--------------------------------
library(pagoo) # Load package
tgz <- system.file('extdata', 'toy_data.tar.gz', package = 'pagoo')
untar(tarfile = tgz, exdir = tempdir()) # Decompress example dataset
files <- list.files(path = tempdir(), full.names = TRUE, pattern = 'tsv$|fasta$') # List files
files

## -----------------------------------------------------------------------------
data_file <- grep("case_df.tsv", files, value = TRUE)
data <- read.table(data_file, header = TRUE, sep = '\t', quote = '')
head(data)

## ---- eval=FALSE--------------------------------------------------------------
#  pg <- pagoo(data = data)

## -----------------------------------------------------------------------------
# Organism metadata
orgs_file <- grep("case_orgs_meta.tsv", files, value = TRUE)
orgs_meta <- read.table(orgs_file, header = TRUE, sep = '\t', quote = '')
head(orgs_meta)

## -----------------------------------------------------------------------------
# Cluster metadata
clust_file <- grep("case_clusters_meta.tsv", files, value = TRUE)
clust_meta <- read.table(clust_file, header = TRUE, sep = '\t', quote = '')
head(clust_meta)

## -----------------------------------------------------------------------------
# Sequences
fasta_files <- grep("[.]fasta", files, value = TRUE) # List fasta files
names(fasta_files) <- sub('[.]fasta', '', basename(fasta_files)) # Name them 
# Read fasta files with Biostrings:
library(Biostrings)
sq <- lapply(fasta_files, readDNAStringSet)
class(sq) # Is list?
length(sq) # One list element per organism
names(sq) # Names are the same as in data$org
class(sq[[1]]) # Class of each element of the list

## ---- message=FALSE-----------------------------------------------------------
p <- pagoo(data = data, # Required data
           org_meta = orgs_meta, # Organism's metadata
           cluster_meta = clust_meta, # Cluster's metadata
           sequences = sq) # Sequences

## -----------------------------------------------------------------------------
host_df <- data.frame(org = p$organisms$org, host = c("Cow", "Dog", "Cat", "Cow", "Sheep"))
p$add_metadata(map = "org", host_df)
p$organisms

## -----------------------------------------------------------------------------
tmp <- paste(tempdir(), "pangenome", sep = "/")
p$write_pangenome(dir = tmp)
list.files(tmp, full.names = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
# Remove unused directory
unlink("pangenome", recursive = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  rds <- paste(tempdir(), "pangenome.RDS", sep = "/")
#  p$save_pangenomeRDS(file = rds)
#  p2 <- load_pangenomeRDS(rds)

