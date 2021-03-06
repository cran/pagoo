#' @name pagoo
#' @title Create a Pagoo Object
#' @description This is the main function to load a pagoo object. It's safer and
#' more friendly than using pagoo's class constructors (\link{PgR6},
#' \link{PgR6M}, and \link{PgR6MS}). This function returns either a
#' \code{\link{PgR6M}} class object, or a \code{\link{PgR6MS}} class object,
#' depending on the parameters provided. If sequences are provided, it returns
#' the latter. See below for more details.
#' @param data A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing at least the
#'   following columns: \code{gene} (gene name), \code{org} (organism name to
#'   which the gene belongs to), and \code{cluster} (group of orthologous to
#'   which the gene belongs to). More columns can be added as metadata for each
#'   gene.
#' @param org_meta (optional) A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing additional
#'   metadata for organisms. This \code{data.frame} must have a column named
#'   "org" with valid organisms names (that is, they should match with those
#'   provided in \code{data}, column \code{org}), and additional columns will be
#'   used as metadata. Each row should correspond to each organism.
#' @param cluster_meta (optional) A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing additional
#'   metadata for clusters. This \code{data.frame} must have a column named
#'   "cluster" with valid organisms names (that is, they should match with those
#'   provided in \code{data}, column \code{cluster}), and additional columns
#'   will be used as metadata. Each row should correspond to each cluster.
#' @param sequences (optional) Can accept: 1) a named \code{list} of named
#'   \code{character} vector. Name of list are names of organisms, names of
#'   character vector are gene names; or 2) a named \code{list} of
#'   \code{\link[Biostrings:XStringSet-class]{DNAStringSetList}} objects (same
#'   requirements as (1), but with BStringSet names as gene names); or 3) a
#'   \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}} (same
#'   requirements as (2) but \code{DNAStringSetList} names are organisms names).
#'   If this parameter is used, then a \code{\link{PgR6MS}} class object is
#'   returned.
#' @param core_level The initial core_level (that's the percentage of organisms a
#'   core cluster must be in to be considered as part of the core genome). Must
#'   be a number between 100 and 85, (default: 95). You can change it later by
#'   using the \code{$core_level} field once the object was created.
#' @param sep A separator. By default is '__'(two underscores). It will be used
#'   to create a unique \code{gid} (gene identifier) for each gene. \code{gid}s
#'   are created by pasting \code{org} to \code{gene}, separated by \code{sep}.
#' @param verbose \code{logical}. Whether to display progress messages when
#'   loading class.
#' @details This package uses
#' [R6](https://r6.r-lib.org/articles/Introduction.html) classes to provide a
#' unified, comprehensive, standardized, but at the same time flexible, way to
#' analyze a pangenome. The idea is to have a single object which contains both
#' the data and the basic methods to analyze them, as well as manipulate fields,
#' explore, and to use in harmony with the already existing and extensive list of R
#' packages available created for comparative genomics and genetics.
#'
#' For more information, tutorials, and resources, please visit https://iferres.github.io/pagoo/ .
#'
#' @section Index:
#' \subsection{Active Bindings}{
#'    \itemize{
#'    \item \href{#field-pan_matrix}{\code{$pan_matrix}}
#'    \item \href{#field-organisms}{\code{$organisms}}
#'    \item \href{#field-clusters}{\code{$clusters}}
#'    \item \href{#field-genes}{\code{$genes}}
#'    \item \href{#field-sequences}{\code{$sequences}}
#'    \item \href{#field-core_level}{\code{$core_level}}
#'    \item \href{#field-core_genes}{\code{$core_genes}}
#'    \item \href{#field-core_clusters}{\code{$core_clusters}}
#'    \item \href{#field-core_sequences}{\code{$core_sequences}}
#'    \item \href{#field-shell_genes}{\code{$shell_genes}}
#'    \item \href{#field-shell_clusters}{\code{$shell_clusters}}
#'    \item \href{#field-shell_sequences}{\code{$shell_sequences}}
#'    \item \href{#field-cloud_genes}{\code{$cloud_genes}}
#'    \item \href{#field-cloud_clusters}{\code{$cloud_clusters}}
#'    \item \href{#field-cloud_sequences}{\code{$cloud_sequences}}
#'    }
#' }
#'
#' @section Active bindings:
#' \describe{\itemize{
#'     \if{html}{\out{<a id="field-pan_matrix"></a>}}
#'     \item{\bold{\code{$pan_matrix}}}{ The panmatrix. Rows are organisms, and
#'     columns are groups of orthologous. Cells indicates the presence (>=1) or
#'     absence (0) of a given gene, in a given organism. Cells can have values
#'     greater than 1 if contain in-paralogs.}
#'     \if{html}{\out{<a id="field-organisms"></a>}}
#'     \item{\bold{\code{$organisms}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with available
#'     organism names, and organism number identifier as \code{rownames()}. (Dropped
#'     organisms will not be displayed in this field, see \code{$dropped} below).
#'     Additional metadata will be shown if provided, as additional columns.}
#'     \if{html}{\out{<a id="field-clusters"></a>}}
#'     \item{\bold{\code{$clusters}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with the groups
#'     of orthologous (clusters). Additional metadata will be shown as additional columns,
#'     if provided before. Each row corresponds to each cluster.}
#'     \if{html}{\out{<a id="field-genes"></a>}}
#'     \item{\bold{\code{$genes}}}{ A \code{\link[IRanges:DataFrameList-class]{SplitDataFrameList}} object with
#'     one entry per cluster. Each element contains a \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'     with gene ids (\code{<gid>}) and additional metadata, if provided. \code{gid} are
#'     created by \code{paste}ing organism and gene names, so duplication in gene names
#'     are avoided.}
#'     \if{html}{\out{<a id="field-sequences"></a>}}
#'     \item{\bold{\code{$sequences}}}{ A \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}} with the
#'     set of sequences grouped by cluster. Each group is accessible as were a list. All
#'     \code{Biostrings} methods are available.}
#'     \if{html}{\out{<a id="field-core_level"></a>}}
#'     \item{\bold{\code{$core_level}}}{ The percentage of organisms a gene must be in
#'     to be considered as part of the coregenome. \code{core_level = 95} by default.
#'     Can't be set above 100, and below 85 raises a warning.}
#'     \if{html}{\out{<a id="field-core_genes"></a>}}
#'     \item{\bold{\code{$core_genes}}}{ Like \code{genes}, but only showing core genes.}
#'     \if{html}{\out{<a id="field-core_clusters"></a>}}
#'     \item{\bold{\code{$core_clusters}}}{ Like \code{$clusters}, but only showing core
#'     clusters.}
#'     \if{html}{\out{<a id="field-core_sequences"></a>}}
#'     \item{\bold{\code{$core_sequences}}}{ Like \code{$sequences}, but only showing core
#'     sequences.}
#'     \if{html}{\out{<a id="field-cloud_genes"></a>}}
#'     \item{\bold{\code{$cloud_genes}}}{ Like \code{genes}, but only showing cloud genes.
#'     These are defined as those clusters which contain a single gene (singletons), plus
#'     those which have more than one but its organisms are probably clonal due to identical
#'     general gene content. Colloquially defined as strain-specific genes.}
#'     \if{html}{\out{<a id="field-cloud_clusters"></a>}}
#'     \item{\bold{\code{$cloud_clusters}}}{ Like \code{$clusters}, but only showing cloud
#'     clusters as defined above.}
#'     \if{html}{\out{<a id="field-cloud_sequences"></a>}}
#'     \item{\bold{\code{$cloud_sequences}}}{ Like \code{$sequences}, but only showing cloud
#'     sequences as defined above.}
#'     \if{html}{\out{<a id="field-shell_genes"></a>}}
#'     \item{\bold{\code{$shell_genes}}}{ Like \code{genes}, but only showing shell genes.
#'     These are defined as those clusters than don't belong neither to the core genome,
#'     nor to cloud genome. Colloquially defined as genes that are present in some but not
#'     all strains, and that aren't strain-specific.}
#'     \if{html}{\out{<a id="field-shell_clusters"></a>}}
#'     \item{\bold{\code{$shell_clusters}}}{ Like \code{$clusters}, but only showing shell
#'     clusters, as defined above.}
#'     \if{html}{\out{<a id="field-shell_sequences"></a>}}
#'     \item{\bold{\code{$shell_sequences}}}{ Like \code{$sequences}, but only showing shell
#'     sequences, as defined above.}
#'     \if{html}{\out{<a id="field-summary_stats"></a>}}
#'     \item{\bold{\code{$summary_stats}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with
#'     information about the number of core, shell, and cloud clusters, as well as the
#'     total number of clusters.}
#'     \if{html}{\out{<a id="field-random_seed"></a>}}
#'     \item{\bold{\code{$random_seed}}}{ The last \code{.Random.seed}. Used for
#'     reproducibility purposes only.}
#'     \if{html}{\out{<a id="field-dropped"></a>}}
#'     \item{\bold{\code{$dropped}}}{ A \code{character} vector with dropped organism
#'     names, and organism number identifier as \code{names()}}
#'  }
#' }
#' @section Methods:
#' \describe{
#'     Below is a comprehensive description of all the methods provided by the object.
#'
#'     \if{html}{\out{<a id="method-add_metadata"></a>}}
#'     \item{\bold{Add metadata}}{
#'             \subsection{Description:}{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$add_metadata(map = 'org', df)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{map}}: \code{character} identifying the metadata to map. Can
#'                     be one of \code{"org"}, \code{"group"}, or \code{"gid"}.}
#'                     \item{\bold{\code{df}}: \code{data.frame} or \code{DataFrame} with the metadata to
#'                     add. For each case, a column named as \code{"map"} must exists, which should
#'                     contain identifiers for each element. In the case of adding gene (\code{gid})
#'                     metadata,each gene should be referenced by the name of the organism and the name
#'                     of the gene as provided in the \code{"data"} data.frame, separated by the
#'                     \code{"sep"} argument.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     \code{self} invisibly, but with additional metadata.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-drop"></a>}}
#'     \item{\bold{Drop an organism}}{
#'             \subsection{Description:}{
#'                    Drop an organism from the dataset. This method allows to hide an organism from
#'                    the real dataset, ignoring it in downstream analyses. All the fields and
#'                    methods will behave as it doesn't exist. For instance, if you decide to drop
#'                    organism 1, the \code{$pan_matrix} field (see below) would not show it when
#'                    called.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$drop(x)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{x}}:  \code{character} or \code{numeric}. The name of the
#'                     organism wanted to be dropped, or its numeric id as returned in
#'                     \code{$organism} field (see below).}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     \code{self} invisibly, but with \code{x} dropped. It isn't necessary
#'                     to assign the function call to a new object, nor to re-write it as R6 objects
#'                     are mutable.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-recover"></a>}}
#'     \item{\bold{Recover a dropped organism}}{
#'             \subsection{Description:}{
#'                   Recover a previously \code{$drop()}ped organism (see above). All fields
#'                   and methods will start to behave considering this organism again.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$recover(x)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{x}}: \code{character} or \code{numeric}. The name of the
#'                     organism wanted to be recover, or its numeric id as returned in
#'                     \code{$dropped} field (see below).}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     \code{self} invisibly, but with \code{x} recovered. It isn't necessary
#'                     to assign the function call to a new object, nor to re-write it as R6 objects
#'                     are mutable.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-write_pangenome"></a>}}
#'     \item{\bold{Write a pangenome as flat (text) files.}}{
#'             \subsection{Description:}{
#'                   Write the pangenome data as flat tables (text). Is not the most recommended way
#'                   to save a pangenome, since you can loose information as numeric precision,
#'                   column classes (factor, numeric, integer), and the state of the object itself
#'                   (i.e. dropped organisms, or core_level), loosing reproducibility. Use
#'                   \code{save_pangenomeRDS} for a more precise way of saving a pagoo object.
#'                   Still, it is useful if you want to work with the data outside R, just keep
#'                   the above in mind.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$write_pangenome(dir = "pangenome", force = FALSE)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{dir}}: The unexisting directory name where to put the data files. Default
#'                     is "pangenome".}
#'                     \item{\bold{\code{force}}: \code{logical}. Whether to overwrite the directory if it already
#'                     exists. Default: \code{FALSE}.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A directory with at least 3 files. "data.tsv" contain the basic
#'                     pangenome data as it is provided to the \code{data} argument in the
#'                     initialization method (\code{$new(...)}). "clusters.tsv" contain any metadata
#'                     associated to the clusters. "organisms.tsv" contain any metadata associated to
#'                     the organisms. The latter 2 files will contain a single column if no metadata
#'                     was provided.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-save_pangenomeRDS"></a>}}
#'     \item{\bold{Save a pangenome as a RDS (binary) file.}}{
#'             \subsection{Description:}{
#'                  Save a pagoo pangenome object. This function provides a method for saving a pagoo
#'                  object and its state into a "RDS" file. To load the pangenome, use the
#'                  \code{load_pangenomeRDS} function in this package. It *should* be compatible between
#'                  pagoo versions, so you could update pagoo and still recover the same pangenome. Even
#'                  \code{sep} and \code{core_level} are restored unless the user provides those
#'                  arguments in \code{load_pangenomeRDS}. \code{dropped} organisms also kept hidden, as
#'                  you where working with the original object.
#'             }
#'
#'             \subsection{Usage:}{
#'                  \verb{                  }\code{$save_pangenomeRDS(file = "pangenome.rds")}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{file}}: The name of the file to save. Default: "pangenome.rds".}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     Writes a list with all the information needed to restore the object by
#'                     using the load_pangenomeRDS function, into an RDS (binary) file.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-clone"></a>}}
#'     \item{\bold{Clone a pagoo object.}}{
#'             \subsection{Description:}{
#'                   The objects of this class are clonable with this method.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$clone(deep = FALSE)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{deep}}: \code{character} identifying the metadata to map. Can
#'                     be one of \code{"org"}, \code{"group"}, or \code{"gid"}.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     Whether to make a deep clone.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-dist"></a>}}
#'     \item{\bold{Compute distances}}{
#'             \subsection{Description:}{
#'                   Compute distance between all pairs of genomes. The default dist method is
#'                   \code{"bray"} (Bray-Curtis distance). Another used distance method is \code{"jaccard"},
#'                   but you should set \code{binary = FALSE} (see below) to obtain a meaningful result.
#'                   See \code{\link[vegan]{vegdist}} for details, this is just a wrapper function.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$dist(
#'                   method = "bray",
#'                   binary = FALSE,
#'                   diag = FALSE,
#'                   upper = FALSE,
#'                   na.rm = FALSE,
#'                   ...
#'                   )}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{method}}: The distance method to use. See \link[vegan]{vegdist}
#'                     for available methods, and details for each one.}
#'                     \item{\bold{\code{binary}}: Transform abundance matrix into a presence/absence
#'                     matrix before computing distance.}
#'                     \item{\bold{\code{diag}}: Compute diagonals.}
#'                     \item{\bold{\code{upper}}: Return only the upper diagonal.}
#'                     \item{\bold{\code{na.rm}}: Pairwise deletion of missing observations when
#'                     computing dissimilarities.}
#'                     \item{\bold{\code{...}}: Other parameters. See \link[vegan]{vegdist} for details.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A \code{dist} object containing all pairwise dissimilarities between genomes.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-pan_pca"></a>}}
#'     \item{\bold{Compute a Principal Component Analysis}}{
#'             \subsection{Description:}{
#'                   Performs a principal components analysis on the panmatrix.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$pan_pca( center = TRUE, scale. = FALSE, ...)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{center}}: a logical value indicating whether the variables should be shifted
#'                     to be zero centered. Alternately, a vector of length equal the number of columns of x can be
#'                     supplied. The value is passed to scale.}
#'                     \item{\bold{\code{scale.}}: a logical value indicating whether the variables should be scaled
#'                     to have unit variance before the analysis takes place. The default is TRUE.}
#'                     \item{\bold{\code{...}}: Other arguments. See \link[stats]{prcomp}}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     Returns a list with class "prcomp". See \link[stats]{prcomp} for more information.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-pg_power_law_fit"></a>}}
#'     \item{\bold{Fit a Power Law Function for the Pangenome}}{
#'             \subsection{Description:}{
#'                   Fits a power law curve for the pangenome rarefaction simulation.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$pg_power_law_fit(raref, ...)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{raref}}: (Optional) A rarefaction matrix, as returned by \code{rarefact()}.}
#'                     \item{\bold{\code{...}}: Further arguments to be passed to \code{rarefact()}. If \code{raref}
#'                     is missing, it will be computed with default arguments, or with the ones provided here.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
#'                     with fitted parameters. An attribute \code{"alpha"} is also returned (If
#'                     \code{alpha>1}, then the pangenome is closed, otherwise is open.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-cg_exp_decay_fit"></a>}}
#'     \item{\bold{Fit an Exponential Decay Function for the Coregenome}}{
#'             \subsection{Description:}{
#'                   Fits an exponential decay curve for the coregenome rarefaction simulation.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$cg_exp_decay_fit(raref, pcounts = 10, ...)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{raref}}: (Optional) A rarefaction matrix, as returned by \code{rarefact()}.}
#'                     \item{\bold{\code{pcounts}}: An integer of pseudo-counts. This is used to better fit the function
#'                     at small numbers, as the linearization method requires to subtract a constant C, which is the
#'                     coregenome size, from \code{y}. As \code{y} becomes closer to the coregenome size, this operation
#'                     tends to 0, and its logarithm goes crazy. By default \code{pcounts=10}.}
#'                     \item{\bold{\code{...}}: Further arguments to be passed to \code{rarefact()}. If \code{raref}
#'                     is missing, it will be computed with default arguments, or with the ones provided here.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
#'                     with fitted intercept and decay parameters.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-fluidity"></a>}}
#'     \item{\bold{Compute Genomic Fluidity}}{
#'             \subsection{Description:}{
#'                   Computes the genomic fluidity, which is a measure of population
#'                   diversity. See \code{\link[micropan]{fluidity}} for more details.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$fluidity(nsim = 10)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{nsim}}:An integer specifying the number of random samples
#'                      to use in the computations.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A list with two elements, the mean fluidity and its sample standard
#'                     deviation over the n.sim computed values.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_barplot"></a>}}
#'     \item{\bold{Plot Accessory Frequency Plot}}{
#'             \subsection{Description:}{
#'                   Plot a barplot with the frequency of genes within the total number of
#'                   genomes.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_barplot()}
#'             }
#'
#'
#'             \subsection{Return:}{
#'                     A barplot, and a \code{gg} object (\code{ggplot2} package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_dist"></a>}}
#'     \item{\bold{Plot a Distance Heatmap}}{
#'             \subsection{Description:}{
#'                   Plot a heatmap showing the computed distance between all pairs of organisms.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_dist(method = "bray", ...)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{method}}: Distance method. One of "Jaccard" (default), or
#'                     "Manhattan",see above.}
#'                     \item{\bold{\code{...}}: More arguments to be passed to \code{\link[micropan]{distManhattan}}}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                      A heatmap (\code{ggplot2::geom_tile()}), and a \code{gg} object (\code{ggplot2}
#'                      package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_binmap"></a>}}
#'     \item{\bold{Plot a Pangenome Binary Map}}{
#'             \subsection{Description:}{
#'                   Plot a pangenome binary map representing the presence/absence of each
#'                   gene within each organism.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_binmap()}
#'             }
#'
#'             \subsection{Return:}{
#'                     A binary map (\code{ggplot2::geom_raster()}), and a \code{gg} object (\code{ggplot2}
#'                     package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_pca"></a>}}
#'     \item{\bold{Plot a PCA}}{
#'             \subsection{Description:}{
#'                   Plot a scatter plot of a Principal Components Analysis.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_pca(colour = NULL, ...))}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{colour}}:The name of the column in \code{$organisms} field
#'                     from which points will take color (if provided). \code{NULL} (default) renders
#'                     black points.}
#'                     \item{\bold{\code{...}}: More arguments to be passed to \code{ggplot2::autoplot()}.}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                     A scatter plot (\code{ggplot2::autoplot()}), and a \code{gg} object
#'                     (\code{ggplot2} package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_pie"></a>}}
#'     \item{\bold{Plot a Pie with Pangenome Categories}}{
#'             \subsection{Description:}{
#'                   Plot a pie chart showing the number of clusters of each pangenome category: core,
#'                   shell, or cloud.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_pie()}
#'             }
#'
#'             \subsection{Return:}{
#'                     A pie chart (\code{ggplot2::geom_bar() + coord_polar()}), and a
#'                     \code{gg} object (\code{ggplot2} package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-gg_curves"></a>}}
#'     \item{\bold{Plot Pangenome Curves}}{
#'             \subsection{Description:}{
#'                   Plot pangenome and/or coregenome curves with the fitted functions returned by
#'                   \code{pg_power_law_fit()} and \code{cg_exp_decay_fit()}. You can add points by
#'                   adding \code{+ geom_points()}, of ggplot2 package.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$gg_curves(what = c("pangenome", "coregenome", ...)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{what}}:  \code{"pangenome"} and/or \code{"coregenome"}.}
#'                     \item{\bold{\code{...}}: ignored}
#'                 }
#'             }
#'             \subsection{Return:}{
#'                    A scatter plot, and a \code{gg} object (\code{ggplot2} package) invisibly.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-runShinyApp"></a>}}
#'     \item{\bold{Run a Shiny App}}{
#'             \subsection{Description:}{
#'                   Launch an interactive shiny app. It contains a sidebar
#'                   with controls and switches to interact with the pagoo
#'                   object. You can drop/recover organisms from the dataset,
#'                   modify the core_level, visualize statistics, plots, and
#'                   browse cluster and gene information. In the main body, it
#'                   contains 2 tabs to switch between summary statistics plots
#'                   and core genome information on one side, and accessory
#'                   genome plots and information on the other.
#'
#'                   The lower part of each tab contains two tables, side by
#'                   side. On the "Summary" tab, the left one contain
#'                   information about core clusters, with one cluster per row.
#'                   When one of them is selected (click), the one on the right
#'                   is updated to show information about its genes (if
#'                   provided), one gene per row. On the "Accessory" tab, a
#'                   similar configuration is shown, but on this case only
#'                   accessory clusters/genes are displayed. There is a slider
#'                   on the sidebar where one can select the accessory
#'                   frequency range to display.
#'
#'                   Give it a try!
#'
#'                   Take into account that big pangenomes can slow down the
#'                   performance of the app. More than 50-70 organisms often
#'                   leads to a delay in the update of the plots/tables.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$runShinyApp()}
#'             }
#'
#'             \subsection{Return:}{
#'                    Opens a shiny app on the browser.
#'             }
#'     }
#'
#'     \if{html}{\out{<a id="method-core_seqs_4_phylo"></a>}}
#'     \item{\bold{Retrieve Core Genes for Phylogeny}}{
#'             \subsection{Description:}{
#'                  A field for obtaining core gene sequences is available (see below),
#'                  but for creating a phylogeny with this sets is useful to: 1) have the
#'                  possibility of extracting just one sequence of each organism on each
#'                  cluster, in case paralogues are present, and 2) filling gaps with
#'                  empty sequences in case the core_level was set below 100\%, allowing
#'                  more genes (some not in 100\% of organisms) to be incorporated to the
#'                  phylogeny. That is the purpose of this special function.
#'             }
#'
#'             \subsection{Usage:}{
#'                   \verb{                  }\code{$core_seqs_4_phylo(max_per_org = 1, fill = TRUE)}
#'             }
#'
#'             \subsection{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{max_per_org}}:  Maximum number of sequences of each organism
#'                     to be taken from each cluster.}
#'                     \item{\bold{\code{fill}}: \code{logical}. If fill \code{DNAStringSet} with
#'                     empty \code{DNAString} in cases where \code{core_level} is set below 100\%,
#'                     and some clusters with missing organisms are also considered.}
#'                 }
#'             }
#'
#'             \subsection{Return:}{
#'                    A \code{DNAStringSetList} with core genes. Order of organisms on each cluster
#'                    is conserved, so it is easier to concatenate them into a super-gene suitable
#'                    for phylogenetic inference.
#'             }
#'     }
#' }
#'
#' @export
pagoo <- function(data, org_meta, cluster_meta, sequences, core_level = 95, sep = "__", verbose = TRUE){

  if (missing(sequences)){

    PgR6M$new(data = data,
              org_meta = org_meta,
              cluster_meta = cluster_meta,
              core_level = core_level,
              sep = sep,
              verbose = verbose)

  } else{

    PgR6MS$new(data = data,
               org_meta = org_meta,
               cluster_meta = cluster_meta,
               sequences = sequences,
               core_level = core_level,
               sep = sep,
               verbose = verbose)
  }
}






#' @name load_pangenomeRDS
#' @title Load A Pagoo Pangenome
#' @description This function loads a pagoo pangenome from a `.RDS` file generated
#' by the `save_pangenomeRDS` method. Objects loaded by this functions keep their
#' states, i.e : dropped/recovered organisms are conserved, as well as the
#' `core_level` setted when the object was originally saved.
#' @param file The path to the pangenome `.RDS` file.
#' @param pkg The package to use to load the object. Shouldn't be necessary to
#' provide, but may be useful in some cases.
#' @param ... Arguments to be passed to the pagoo object. \code{sep} and
#' \code{core_level} overwrite the values stored in the file.
#' @return A \code{PgR6MS} class object, or a \code{PgR6M} object (with or
#' without sequences, respectively).
#' @importFrom utils installed.packages
#' @export
load_pangenomeRDS = function(file, pkg, ...){

  args <- readRDS(file)
  atrs <- attributes(args)
  dots <- list(...)

  if (!"package" %in% names(atrs)) stop("Not recognized rds file.")
  if (atrs$package != "pagoo") stop("Not recognized rds file.")
  if (missing(pkg)){
    pkg <- atrs$parent_package
    if (is.null(pkg)) {
      warning("'$parent_package' attribute not found. Setting 'pagoo' as parent package.", immediate. = TRUE)
      pkg <- "pagoo"
    }
  }
  if (! pkg %in% rownames(installed.packages())) stop(paste(pkg, "not installed."))
  clss <- atrs$class

  dropped <- args$dropped
  args$dropped <- NULL

  sep1 <- args$sep
  sep2 <- dots$sep
  if (!is.null(sep2)){
    if (sep2 != sep1) warning('"sep" argument provided is not the same as the original object. Overwriting.')
    args$sep <- sep2
  }

  core_level1 <- args$core_level
  core_level2 <- dots$core_level
  if (!is.null(core_level2)){
    if (core_level2 != core_level1) warning('"core_level" argument provided is not the same as the original object. Overwriting.')
    args$core_level <- core_level2
  }


  i <- 1L
  clp <- NA_character_
  while (! clp[1] %in% clss ) {
    cmd <- paste0(pkg, ":::", clss[i], "$new")
    p <- try(do.call(eval(parse(text = cmd)), args), silent = TRUE)
    clp <- class(p)
    if ( ! clp[1] %in% clss ){
      wrn <- paste(clss[i], "not in", pkg, "namespace. Rolling back to", clss[i+1],".")
      warning(wrn, immediate. = TRUE)
    }
    i <- i + 1L
  }

  if (!is.null(dropped)){
    message(paste("Dropping the following organisms:", paste0(dropped, collapse = " "), collapse = " "))
    p$drop(dropped)
  }

  return(p)
}
