#Necessary for S4 to recognize "data.table" as a class
setOldClass(c("data.frame"))
setOldClass(c("data.table", "data.frame"))

#' The primary object for BAGEL that contains all samples and results
#'
#' @slot samples A list of sample objects containing sample-level information
#' @slot prop_table Summary tables with motif counts normalized by sample sums
#' @slot counts_table Summary tables with unnormalized motif counts
#' @export
setClass("bagel", representation(samples = "data.table", prop_table = "matrix",
                                 counts_table = "matrix"),
         prototype(samples = NULL, prop_table = matrix(nrow = 0, ncol = 0),
                   counts_table = matrix(nrow = 0, ncol = 0)))

setMethod("show", "bagel",
          function(object)cat("BayeSig Object containing \nSamples: ",
                              if (!is.null(object@samples)) {
                                noquote(paste(length(unique(
                                  object@samples$Tumor_Sample_Barcode)),
                                  collapse = ", "))
                                }else{
                                  "Empty"
                                    },
                              "\nCounts Table Rows: ",
                              if (!is.null(object@counts_table)) {
                                rowSums(object@counts_table)
                                }else{
                                  "Empty"
                                  })
)

#' Set samples for bagel object
#'
#' @param bay Bagel object we input sample into
#' @param samp Sample DataFrame
#' @return Sets sample slot {no return}
#' @examples
#' samp <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' set_samples(bay, samp)
#' @export
set_samples <- function(bay, samp) {
  eval.parent(substitute(bay@samples <- samp))
}

#' Return samples names for bagel object
#'
#' @param bay Bagel object containing samples
#' @return Returns names of samples in bagel object
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' sample_names(bay)
#' @export
sample_names <- function(bay) {
  return(unique(bay@samples$Tumor_Sample_Barcode))
}

#' Return sample from bagel object
#'
#' @param bay Bagel object containing samples
#' @param sample_name Sample name to subset by
#' @return Returns sample dataframe subset to a single sample
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' subset_samples(bay, "public_LUAD_TCGA-97-7938.vcf")
#' @export
subset_samples <- function(bay, sample_name) {
  return(bay@samples[which(bay@samples$Tumor_Sample_Barcode == sample_name), ])
}

# Result object/methods -------------------------------

#' Object containing deconvolved/predicted signatures, sample weights, and
#' the bagel object the result was generated from
#'
#' @slot signatures A matrix of signatures by mutational motifs
#' @slot samples A matrix of samples by signature weights
#' @slot type Describes how the signatures/weights were generated
#' @slot bagel The bagel object the results were generated from
#' @export
setClass("Result", representation(signatures = "matrix", samples = "matrix",
                                  type = "character", bagel = "bagel"))

#' Return sample from bagel object
#'
#' @param result Result object containing signatures and weights
#' @param name_vector Vector of user-defined signature names
#' @return Result object with user-defined signatures names
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' name_signatures(result, c("smoking", "uv", "apobec", "unknown"))
#' @export
name_signatures <- function(result, name_vector) {
  num_sigs <- length(colnames(result@signatures))
  if (length(name_vector) != num_sigs) {
    stop(paste("Please provide a full list of signatures names (length = ",
               num_sigs, ")", sep = ""))
  }
  eval.parent(substitute(colnames(result@signatures) <- name_vector))
}
