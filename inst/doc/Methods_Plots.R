## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(pagoo) # Load package
toy_rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(toy_rds)

## -----------------------------------------------------------------------------
p$dist()

## ---- eval=FALSE--------------------------------------------------------------
#  p$dist(method = "jaccard", binary = TRUE)

## -----------------------------------------------------------------------------
pca <- p$pan_pca()
summary(pca)

## ---- message=FALSE-----------------------------------------------------------
library(ggplot2)
library(patchwork) # To arrange plots

## -----------------------------------------------------------------------------
# Basic
pie1 <- p$gg_pie() + ggtitle("Default")

# Customize with ggplot2
pie2 <- pie1 + 
  ggtitle("Customized") + 
  theme_bw(base_size = 15) + 
  scale_fill_brewer(palette = "Blues")

# Arrange (patchwork) and plot
pie1 + pie2

## -----------------------------------------------------------------------------
p$gg_barplot()

## -----------------------------------------------------------------------------
p$gg_binmap()

## -----------------------------------------------------------------------------
p$gg_curves()

## ---- warning=FALSE-----------------------------------------------------------
p$gg_curves() + 
  ggtitle("Pangenome and Coregenome curves") + 
  geom_point() + 
  facet_wrap(~Category, scales = 'free_y') + 
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Accent")

## -----------------------------------------------------------------------------
p$organisms

## ---- warning=FALSE-----------------------------------------------------------
p$gg_pca(colour = 'host', size = 4) + 
  theme_bw(base_size = 15) +
  scale_color_brewer(palette = "Set2")

## ---- eval=FALSE--------------------------------------------------------------
#  p$runShinyApp()

