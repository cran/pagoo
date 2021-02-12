## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(pagoo) # Load package
toy_rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(toy_rds)

## -----------------------------------------------------------------------------
p$clusters
p$core_clusters
p$shell_clusters
p$cloud_clusters

## -----------------------------------------------------------------------------
# From clusters 1 to 3
p[c(1, 5, 18)]$pan_matrix
p[c(1, 5, 18)]$genes
p[c(1, 5, 18)]$clusters
p[c(1, 5, 18)]$sequences
p[c(1, 5, 18)]$organisms

## -----------------------------------------------------------------------------
shell_clust <- p$shell_clusters$cluster[1] # [1] "OG0005"
p[shell_clust]$organisms

## -----------------------------------------------------------------------------
p$pan_matrix[1:3, c(1, 5, 18)]

## -----------------------------------------------------------------------------
p[1:3, c(1, 5, 18)]$pan_matrix # The same as above
p[1:3, c(1, 5, 18)]$genes
p[1:3, c(1, 5, 18)]$clusters
p[1:3, c(1, 5, 18)]$sequences
p[1:3, c(1, 5, 18)]$organisms

## -----------------------------------------------------------------------------
p$organisms
p$summary_stats

## -----------------------------------------------------------------------------
p$drop(1:2) 
p$organisms # Updated organisms !!
p$summary_stats # Updated stats !!

## -----------------------------------------------------------------------------
p$dropped
p$recover( p$dropped )
p$organisms

