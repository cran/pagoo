## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(pagoo) # Load package
toy_rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(toy_rds)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(ggplot2)
library(patchwork)

# 1. Pangenome curves
curves <- p$gg_curves() +                                     # Plot core- and pan-genome curves
          scale_color_manual(values = c('black', 'black')) +  # Customize line colors
          geom_point(alpha = .05, size = 4, color = 'grey') + # Add semi-transparent data points
          theme_bw(base_size = 15) +                          # Customize background theme
          theme(legend.position = 'none',                     # Remove legend
		axis.title = element_text(size = 12),                     # Customize axis title
		axis.text = element_text(size = 12))                      # Customize axis text size

# 2. Gene frequency bar plots
bars <- p$gg_barplot() +                                       # Plot gene frequency distribution
  theme_bw(base_size = 15) +                                   # Customize background color
  theme(axis.title = element_text(size = 12),                  # Customize axis label size
        axis.text=element_text(size = 12)) +                   # Customize axis text size
  geom_bar(stat = 'identity', color = 'black', fill = 'black') # Customize bar color and borders

# 3. PCA of accessory genes colored by host
pca <- p$gg_pca(colour = 'host', size = 4) +                # Plot PCA, color by host
  theme_bw(base_size = 15) +                                # Customize background theme
  theme(legend.position = 'bottom') +                       # Customize legend position
  theme(axis.title = element_text(size = 12),               # Customize axis title
        axis.text = element_text(size = 12))                # Customize axis text size

# 4. Pie chart of core and accessory genes
pie <- p$gg_pie() +                                         # Plot pie chart
  theme_bw(base_size = 15) +                                # Customize background theme
  scale_fill_discrete(guide = guide_legend(keywidth = .75,
                                           keyheight = .75)) + # Customize fill
  scale_fill_brewer(palette = "Blues") +                    # Customize fill color
  scale_x_discrete(breaks = c(0, 25, 50, 75)) +             # Customize axis scales
  theme(legend.position = 'bottom',                         # Customize legend position
        legend.title = element_blank(),                     # Remove legend title
        legend.text = element_text(size = 10),              # Change legend text size
        legend.margin = margin(0, 0, 13, 0),                # Change legend margins
        legend.box.margin = margin(0, 0, 5, 0),             # Change box margins
        axis.title.x = element_blank(),                     # Remove X-axis title
        axis.title.y = element_blank(),                     # Remove Y-axis title
        axis.ticks = element_blank(),                       # Remove axis ticks
        axis.text.x = element_blank())                      # Remove X-axis text


# 5. Use patchwork to arrange plots using math operators
(curves + bars) / (pca + pie)

## ---- eval=FALSE--------------------------------------------------------------
#  library(micropan)
#  
#  micropan::binomixEstimate(p$pan_matrix)

## ---- eval=FALSE--------------------------------------------------------------
#  library(micropan)
#  library(magrittr)
#  
#  p$pan_matrix %>% micropan::fluidity(n.sim = 100)

## ---- eval=FALSE--------------------------------------------------------------
#  # Load required packages
#  if (!require(rBlast)) devtools::install_github('mhahsler/rBLAST')
#  library(Biostrings)
#  library(rBLAST)
#  library(magrittr)
#  
#  db_path <- 'path/to/custom/blastpdb'        # Path to custom blastp db
#  db <- blast(db = db_path, type = 'blastp')  # Set blastdb
#  
#  blast_result <- p$sequences %>%             # Pangenome sequences
#    lapply('[[', 1L) %>%                      # Subset 1 sequence from each cluster
#    Biostrings::DNAStringSet() %>%            # Transform list to DNAStringSet
#    Biostrings::translate() %>%               # Translate DNAStringSet
#    rBLAST::predict.BLAST(db, .)              # Run blastp. Returns data.frame
#  
#  ### ADD METADATA STEP MISSING

## ---- eval=FALSE--------------------------------------------------------------
#  # Load required packages
#  library(magrittr)
#  library(DECIPHER)
#  library(pegas)
#  library(ape)
#  
#  p$core_level <- 100                       # Set core_level to 100% to avoid
#                                            #  AlignTranslation() errors.
#  
#  tajimaD <- p$core_seqs_4_phylo() %>%      # Core genome sequences
#    lapply(DECIPHER::AlignTranslation) %>%  # Align translation
#    lapply(ape::as.DNAbin) %>%              # Transform class to DNAbin
#    lapply(pegas::tajima.test) %>%          # Compute Tajima's test
#    sapply('[[', 'D')                       # Get Tajima's 'D' statistic from each
#  
#  # Which are neutral?
#  which(tajimaD <= 0.2 & tajimaD >= -0.2)

## ---- eval=FALSE--------------------------------------------------------------
#  # Load required packages
#  library(magrittr)
#  library(DECIPHER)
#  library(Biostrings)
#  library(phangorn)
#  library(ggtree)
#  
#  phy <- p$core_seqs_4_phylo() %>%                 # Core genome sequences
#    lapply(DECIPHER::AlignSeqs) %>%                # Align
#    do.call(Biostrings::xscat, .) %>%              # Concatenate alignments
#    setNames(p$organisms$org) %>%                  # Set sequence names
#    as('matrix') %>%                               # Transform to matrix
#    phangorn::phyDat(type = 'DNA') %>%             # Transform to phangorn's phyDat
#    phangorn::dist.ml() %>%                        # Compute distance
#    phangorn::NJ() %T>% {                          # Compute NJ, and assign "phy"
#      {
#        ggtree::ggtree(.) %<+%                        # Create ggtree
#          as.data.frame(p$organisms) +                # Get organisms metadata
#          ggtree::geom_tippoint(aes(colour = host)) + # Add coloured tip points
#          scale_color_brewer(palette = 'Set1')        # Set color palette
#      } %>%
#        print()
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  # Load required packages
#  library(magrittr)
#  library(DECIPHER)
#  library(Biostrings)
#  library(phangorn)
#  library(ggtree)
#  
#  phy <- p$core_seqs_4_phylo() %>%                 # Core genome sequences
#    lapply(DECIPHER::AlignSeqs) %>%                # Align
#    do.call(Biostrings::xscat, .) %>%              # Concatenate alignments
#    setNames(p$organisms$org) %>%                  # Set sequence names
#    as('matrix') %>%                               # Transform to matrix
#    phangorn::phyDat(type = 'DNA') %T>%            # Transform to phangorn's phyDat
#    assign('dat', ., .GlobalEnv) %>%               # Assign to "dat" in .GlobalEnv
#    phangorn::dist.ml() %>%                        # Compute distance
#    phangorn::NJ() %>%                             # Compute NJ (initial tree)
#    phangorn::pml(data = dat, k = 4) %>%           # Compute likelihood with 4 discrete
#                                                   #  gamma distributions.
#    phangorn::optim.pml(rearrangement = "stochastic", # Optimize likelihood with
#                                                   # stochastic rearrangements,
#                        optGamma = TRUE,           # optimize gamma rate parameter,
#                        optInv = TRUE,             # optimize prop of variable size,
#                        model ="GTR") %>%          # and use "GTR" model.
#    magrittr::extract2("tree") %T>% {              # Extract the tree only, and pass
#      {                                            #  it to ggtree.
#        ggtree::ggtree(.) %<+%                        # Create ggtree
#          as.data.frame(p$organisms) +                # Get organisms metadata
#          ggtree::geom_tippoint(aes(colour = host)) + # Add coloured tip points
#          scale_color_brewer(palette = 'Set1')        # Set color palette
#      } %>%
#      print()
#    }

## ---- eval=FALSE--------------------------------------------------------------
#  library(magrittr)
#  library(DECIPHER)
#  library(rhierbaps)
#  library(ape)
#  library(phangorn)
#  
#  # 0. Always use core_level at 100% when using DECIPHER::AlignTranslation()
#  p$core_level <- 100
#  
#  # 1. Align translation of core genes
#  ali <- p$core_seqs_4_phylo() %>%           # Core genome sequences
#    lapply(DECIPHER::AlignTranslation)       # Align translation
#  
#  # 2. Identify neutral core clusters
#  tajD <- ali %>%
#    lapply(ape::as.DNAbin) %>%              # Transform class to DNAbin
#    lapply(pegas::tajima.test) %>%          # Compute Tajima's test
#    sapply('[[', 'D')                       # Subset D statistic
#  neutral <- which(tajD <= 2 & tajD >= -2)
#  
#  # 3. Concatenate neutral core clusters
#  concat_neu <- ali[neutral] %>%            # Select neutral clusters
#    do.call(Biostrings::xscat, .) %>%       # Concatenate alignments
#    setNames(p$organisms$org) %>%           # Set sequence names
#    as('matrix') %>%                        # Transform to matrix
#    tolower()                               # Translate to lower case
#  
#  # 4. Compute structure
#  rhb <- hierBAPS(snp.matrix = concat_neu, # Input matrix alignment
#                  n.pops = 10,             # Max number of subpopulations
#                  max.depth = 1,           # Max depth for hierarchical clustering
#                  n.extra.rounds = 5)      # Extra rounds to ensure convergence
#  
#  # 5. Add lineage as metadata to organisms in pagoo object
#  res <- rhb$partition.df
#  lin <- data.frame(org = as.character(res[, 1]),
#                    lineage = as.factor(res[, 2]))
#  p$add_metadata(map = 'org', data = lin)
#  
#  # 6. Compute tree and plot it with lineage information
#  concat_neu %>%
#    phangorn::phyDat(type = 'DNA') %>%                # Transform to phangorn's phyDat
#    phangorn::dist.ml() %>%                           # Compute distance
#    phangorn::NJ() %>%                                # Compute NJ
#    ggtree::ggtree() %<+%                             # Create ggtree
#       as.data.frame(p$organisms) +                   # Get organisms metadata
#       ggtree::geom_tippoint(aes(colour = lineage))   # Colour tips with lineage info

## ---- eval=FALSE--------------------------------------------------------------
#  library(magrittr)
#  library(IRanges)
#  library(Biostrings)
#  library(DECIPHER)
#  library(ape)
#  library(pegas)
#  
#  # Set core level to 100%. This recipe only works if this is set to 100%.
#  p$core_level <- 100
#  
#  # Create pairs matrix
#  pairs <- data.frame(t(combn(nrow(p$organisms), 2)))
#  colnames(pairs) <- c('org1', 'org2')
#  
#  # Compute paired jaccard similarity, transform to matrix
#  jaccard_sim <- as.matrix(p$dist(method = "jaccard", binary = TRUE))
#  # Fill results matrix
#  pairs$jaccard_sim <- apply(pairs, 1, function(i){
#    ii <- i[1]
#    jj <- i[2]
#    jaccard_sim[ii, jj]
#  })
#  
#  # Return only synonymous polymorphic sites.
#  # First, it removes non-synonymous codons, and then retains only
#  # polymorphyc sites.
#  syn_poly_sites <- p$core_seqs_4_phylo() %>%
#    lapply(function(x){
#      lns <- elementNROWS(x)                              # Align translatation filtering
#      tali <- x[which(lns != 0)] %>%                      #  truncated codons and returning
#        Biostrings::subseq(1L, lns %/% 3 * 3) %>%         #  both DNA and AA alignments.
#        DECIPHER::AlignTranslation(type = "both")
#      syno <-  tali[[2]] %>%                              # Identify non-synonymous
#        Biostrings::consensusMatrix() %>%                 #  codons.
#        magrittr::equals(0) %>%
#        magrittr::not() %>%
#        colSums() %>%
#        magrittr::equals(1) %>%
#        which()
#      neut <- tali[[1]] %>%                               # Remove non-synonymous codons.
#        lapply(function(x){
#          IRanges::successiveViews(
#            x, rep.int(3L, length(x) %/% 3L))
#        }) %>%
#        lapply('[', syno) %>%
#        lapply(unlist) %>%
#        Biostrings::DNAStringSet()
#      poly <- neut %>%                                   # Identify polymorphic sites.
#        Biostrings::consensusMatrix() %>%
#        magrittr::equals(0) %>%
#        magrittr::not() %>%
#        colSums() %>%
#        magrittr::is_greater_than(1) %>%
#        which()
#      lapply(neut, '[', poly) %>%                        # Retain only polymorphic
#        Biostrings::DNAStringSet()                       #  sites
#    }) %>%
#    do.call(Biostrings::xscat, .) %>%                    # Concatenate.
#    setNames(p$organisms$org) %>%                        # Set names.
#    ape::as.DNAbin()                                     # Convert to DNAbin class.
#  
#  # Compute paired nucleotide diversity, and fill results matrix
#  pairs$nuc_div <- apply(pairs, 1, function(i){
#    ii <- i[1]
#    jj <- i[2]
#    pair <- c(syn_poly_sites[ii], syn_poly_sites[jj])
#    pegas::nuc.div(pair)
#  })
#  
#  # Plot correlation with R-base graphics
#  plot(jaccard_sim ~ nuc_div, pairs)
#  abline(lm(jaccard_sim ~ nuc_div, pairs))

