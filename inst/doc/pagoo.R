## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------
library(pagoo, quietly = TRUE, warn.conflicts = FALSE) # Load package
rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
campy <- load_pangenomeRDS(rds) # Load pangenome

## -----------------------------------------------------------------------------
campy$organisms

## -----------------------------------------------------------------------------
campy$clusters

## ---- eval=FALSE--------------------------------------------------------------
#  campy$genes

## -----------------------------------------------------------------------------
campy$dist(method = "bray")

## -----------------------------------------------------------------------------
campy$gg_barplot()

## -----------------------------------------------------------------------------
campy$sequences

