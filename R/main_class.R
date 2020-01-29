#Necessary for S4 to recognize "data.table" as a class
setOldClass(c("data.frame"))
setOldClass(c("data.table", "data.frame"))

#' The primary object for BAGEL that contains all samples and tables
#'
#' @slot variants Data.table of variants and variant-level information
#' @slot counts_table Summary tables with unnormalized motif counts
#' @export
setClass("bagel", representation(variants = "data.table", counts_table =
                                   "matrix", sample_annotations = "data.table"),
         prototype(variants = data.table::data.table(), counts_table = matrix(),
                   sample_annotations = data.table::data.table()))

setMethod("show", "bagel",
          function(object)cat(cat("BAGEL Object containing \n**Variants: \n"),
                              if (!all(is.na(object@variants))) {
                                cat(show(object@variants))
                                }else{
                                  cat("Empty")
                                    },
                              cat("\n**Counts Table Dim and Subset: \n"),
                              if (!all(is.na(object@counts_table))) {
                                cat("Dim: \n")
                                cat(show(dim(object@counts_table)),
                                    "\nSubset Results:\n")
                                cat(show(object@counts_table[seq_len(5),
                                  seq_len(min(3, nrow(object@counts_table))),
                                                                   drop =
                                    FALSE]))
                                }else{
                                  cat("Empty")
                                  },
                              cat("\n**Sample Level Annotations: \n"),
                              if (!all(is.na(object@sample_annotations))) {
                                cat(show(object@sample_annotations))
                              }else{
                                cat("Empty")
                              })
)

# Variant-Level object/methods -------------------------------

#' Set variants table for bagel object
#'
#' @param bay Bagel object we input sample into
#' @param variant_dt Variant DataFrame
#' @return Sets variant slot {no return}
#' @examples
#' variants <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' set_variants(bay, variants)
#' @export
set_variants <- function(bay, variant_dt) {
  eval.parent(substitute(bay@variants <- variant_dt))
}

#' Return sample from bagel object
#'
#' @param bay Bagel object containing samples
#' @param sample_name Sample name to subset by
#' @return Returns sample dataframe subset to a single sample
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' subset_variants_by_samples(bay, "public_LUAD_TCGA-97-7938.vcf")
#' @export
subset_variants_by_samples <- function(bay, sample_name) {
  return(bay@variants[which(bay@variants$Tumor_Sample_Barcode == sample_name),
                      ])
}

# Sample-Level object/methods -------------------------------

#' Set sample level annotations for bagel object
#'
#' @param bay Bagel object we input sample into
#' @param annotations Sample DataFrame
#' @return Sets sample_annotations slot {no return}
#' @examples
#' annotations <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' set_sample_annotations(bay, annotations)
#' @export
set_sample_annotations <- function(bay, annotations) {
  eval.parent(substitute(bay@sample_annotations <- annotations))
}

#' Initialize sample annotation data.table with sample names from variants
#'
#' @param bay Bagel object we input sample into
#' @return Sets sample_annotations slot {no return}
#' @examples
#' annotations <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' init_sample_annotations(bay)
#' @export
init_sample_annotations <- function(bay) {
  #samples <- unique(tools::file_path_sans_ext(
  #  bay@variants$Tumor_Sample_Barcode))
  samples <- unique(bay@variants$Tumor_Sample_Barcode)
  sample_dt <- data.table::data.table(Samples = samples)
  eval.parent(substitute(bay@sample_annotations <- sample_dt))
}

#' Adds sample annotation to bagel object with available samples
#'
#' @param bay Bagel object we input sample into
#' @return Sets sample_annotations slot {no return}
#' @examples
#' annotations <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' init_sample_annotations(bay)
#' @export
add_sample_annotations <- function(bay, annotations, sample_column = Sample_ID,
                                   columns_to_add) {
  bay_annotations <- get_sample_annotations(bay)
  if (all(is.na(bay_annotations))) {
    stop(strwrap(prefix = " ", initial = "", "Please run init_sample_annotations
                 on this bagel before adding sample annotations."))
  }
  if (!sample_column %in% colnames(annotations)) {
    stop(strwrap(prefix = " ", initial = "", "User-defined sample_column is
                 not in input annotations, please check and rerun."))
  }
  if (!all(bay_annotations$Samples %in%
          annotations[, sample_column])) {
    stop(strwrap(prefix = " ", initial = "", "Some samples are missing
                 annotations, please check input annotations and rerun."))
  }
  if (!all(columns_to_add %in% colnames(annotations))) {
    stop(strwrap(prefix = " ", initial = "", paste("Some user-defined
                                                   columns_to_add are not in
                                                   the input annotations, (",
                 toString(columns_to_add[which(!columns_to_add %in%
                                        colnames(annotations))]),
                 ") please check and rerun.", sep = "")))
  }
  matches <- match(bay_annotations$Samples,
                  annotations[, sample_column])
  bay_annotations <- cbind(bay_annotations, annotations[matches, columns_to_add,
                                                       drop = FALSE])
  eval.parent(substitute(bay@sample_annotations <- bay_annotations))
}

#' Return sample annotation from bagel object
#'
#' @param bay Bagel object we input sample into
#' @return Sets sample_annotations slot {no return}
#' @examples
#' annotations <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' init_sample_annotations(bay)
#' show_sample_annotations(bay)
#' @export
get_sample_annotations <- function(bay) {
  return(bay@sample_annotations)
}

#' Return samples names for bagel object
#'
#' @param bay Bagel object containing samples
#' @return Returns names of samples in bagel object
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' sample_names(bay)
#' @export
get_sample_names <- function(bay) {
  return(unique(bay@variants$Tumor_Sample_Barcode))
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
                                  type = "character", bagel = "bagel",
                                  log_lik = "numeric"))

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
  eval.parent(substitute(rownames(result@samples) <- name_vector))
}

# Result Grid object/methods -------------------------------

#' Object containing the result objects generated from the combination of
#' annotations and a range of k values
#'
#' @slot grid_params The parameters the result grid was created using
#' @slot result_list A list of result objects with different parameters
#' @slot grid_table A summary table of the result objects in result_list
#' @export
setClass("Result_Grid", representation(grid_params = "data.table",
                                       result_list = "list",
                                       grid_table = "data.table"))
