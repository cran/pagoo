## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE,  warning=FALSE, message=FALSE--------------------------------
library(pagoo) # Load package
rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(rds)

## -----------------------------------------------------------------------------
p$summary_stats

## -----------------------------------------------------------------------------
p$core_level       
p$core_level <- 100 # Change value
p$summary_stats     # Updated object

## -----------------------------------------------------------------------------
p$core_level <- 95

## -----------------------------------------------------------------------------
p$pan_matrix[, 1:5]

## -----------------------------------------------------------------------------
p$genes

## ---- eval=FALSE--------------------------------------------------------------
#  unlist(p$genes, use.names = FALSE)

## -----------------------------------------------------------------------------
p$clusters

## -----------------------------------------------------------------------------
p$sequences             # List all sequences grouped by cluster
p$sequences[["group0001"]]  # List first cluster

## -----------------------------------------------------------------------------
p$organisms

