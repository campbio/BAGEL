#Necessary for S4 to recognize "data.table" as a class
setOldClass(c("data.frame"))
setOldClass(c("data.table", "data.frame"))

# Count Tables object/methods -------------------------------

#' Object containing the result objects generated from the combination of
#' annotations and a range of k values
#'
#' @slot table_list A list of data.tables with counts
#' @slot table_name The parameters the result grid was created using
#' @slot description A summary table of the result objects in result_list
#' a list of lists. The nested lists created combined (rbind) tables, and the
#' tables at the first list level are modelled independantly. Combined tables
#' must be named.
#' list("tableA", comboTable = list("tableC", "tableD"))
#' @export
setClass("Count_Tables", representation(table_list = "list",
                                        table_name = "list",
                                        description = "list"))

create_count_table <- function(bay, table, name, description = NA,
                               return_instead = FALSE) {
  tab <- bay@count_tables

  #Check that table names are unique
  if (name %in% tab@table_name) {
    stop(paste("Table names must be unique. Current table names are: ",
               paste(tab@table_name, collapse = ", "), sep = ""))
  }

  if (is(table, "matrix")) {
    tab@table_list[[name]] <- table
  } else {
    stop("Please provide only a table or a compounding for tables")
  }

  tab@table_name[[name]] <- name
  tab@description[[name]] <- description
  if (return_instead) {
    return(tab)
  } else {
    eval.parent(substitute(bay@count_tables <- tab))
  }
}

create_variant_table <- function(bay, variant_annotation, name,
                                 description = NA, data_factor = NA,
                                 return_instead = FALSE) {
  tab <- bay@count_tables
  variants <- bay@variants

  #Check that table names are unique
  if (name %in% tab@table_name) {
    stop(paste("Table names must be unique. Current table names are: ",
               paste(tab@table_name, collapse = ", "), sep = ""))
  }

  #Check that variant column exists
  if (variant_annotation %in% colnames(variants)) {
    column_data <- variants[[variant_annotation]]
    sample_names <- unique(variants$Tumor_Sample_Barcode)
    num_samples <- length(sample_names)
    default_factor <- levels(factor(column_data))
    variant_tables <- vector("list", length = num_samples)
    for (i in seq_len(num_samples)) {
      sample_index <- which(variants$Tumor_Sample_Barcode == sample_names[i])
      if (!is.na(data_factor)) {
        variant_tables[[i]] <- table(factor(column_data[sample_index],
                                      levels = data_factor))
      } else {
        variant_tables[[i]] <- table(factor(column_data[sample_index],
                                      levels = default_factor))
      }
    }
    table <- do.call(cbind, variant_tables)
    colnames(table) <- sample_names
  } else {
    stop(paste("That variant annotation does not exist,",
               " existing annotations are: ", paste(colnames(variants),
                                                    collapse = ", "), sep = ""))
  }

  tab@table_list[[name]] <- table
  tab@table_name[[name]] <- name
  tab@description[[name]] <- description
  if (return_instead) {
    return(tab)
  } else {
    eval.parent(substitute(bay@count_tables <- tab))
  }
}

combine_count_tables <- function(bay, to_comb, name, description = NA) {
  tab <- bay@count_tables

  #Check that table names are unique
  if (name %in% tab@table_name) {
    stop(paste("Table names must be unique. Current table names are: ",
               paste(tab@table_name, collapse = ", "), sep = ""))
  }

  if (all(to_comb %in% tab@table_name)) {
    combo_table <- NULL
    for(i in 1:length(to_comb)) {
      combo_table <- rbind(combo_table, tab@table_list[[to_comb[i]]])
    }
    tab@table_list[[name]] <- combo_table
  } else {
    stop(paste("User specified table: ",
               setdiff(to_comb, tab@table_name), " does not exist, please ",
                       "create prior to creating compound table. ",
               "Current table names are: ", paste(tab@table_name,
                                                  collapse = ", "), sep = ""))
  }
  tab@table_name[[name]] <- name
  tab@description[[name]] <- description
  eval.parent(substitute(bay@count_tables <- tab))
}

semi_split <- function(compound) {
  return(unlist(strsplit(compound, ";")))
}

drop_count_table <- function(bay, table_name) {
  tab <- bay@count_tables
  if (!table_name %in% tab@table_name) {
    stop(paste(table_name, " does not exist. Current table names are: ",
               tab@table_name, sep = ""))
  }
  tab@table_list[[table_name]] <- NULL
  tab@table_name[[table_name]] <- NULL
  tab@description[[table_name]] <- NULL
  eval.parent(substitute(bay@count_tables <- tab))
}

setMethod("show", "Count_Tables",
          function(object)cat("Count_Tables Object containing: ",
                              "\n**Count Tables: \n",
                              apply(cbind(do.call("rbind", lapply(
                                object@table_list, dim)), "\n"), 1, paste),
                              "\n**Names: \n",
                              paste(unlist(object@table_name), "\n", sep = ""),
                              "\n**Descriptions: \n",
                              paste(unlist(object@description), "\n", sep = ""))
)

# Primary bagel object/methods -------------------------------

#' The primary object for BAGEL that contains all samples and tables
#'
#' @slot variants Data.table of variants and variant-level information
#' @slot count_tables Summary table with per-sample unnormalized motif counts
#' @slot sample_annotations Sample-level annotations (e.g. age, sex, primary)
#' @export
setClass("bagel", representation(variants = "data.table", count_tables =
                                   "Count_Tables",
                                 sample_annotations = "data.table"),
         prototype(variants = data.table::data.table(),
                   count_tables = new("Count_Tables"),
                   sample_annotations = data.table::data.table()))

setMethod("show", "bagel",
          function(object)cat(cat("BAGEL Object containing \n**Variants: \n"),
                              if (!all(is.na(object@variants))) {
                                cat(methods::show(object@variants))
                                }else{
                                  cat("Empty")
                                    },
                              cat("\n**Count_Tables Object containing: \n"),
                              if (length(object@count_tables@table_name) > 0) {
                                cat("\n**Count Tables: \n",
                                    apply(cbind(do.call("rbind", lapply(
                                      object@count_tables@table_list, dim)),
                                      "\n"), 1, paste),
                                    "\n**Names: \n", paste(
                                        unlist(object@count_tables@table_name),
                                        "\n", sep = ""), "\n**Descriptions: \n",
                                    paste(unlist(
                                      object@count_tables@description), "\n",
                                      sep = ""))
                                }else{
                                  cat("Empty")
                                  },
                              cat("\n**Sample Level Annotations: \n"),
                              if (!all(is.na(object@sample_annotations))) {
                                cat(methods::show(object@sample_annotations))
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
#' @param annotations table of sample-level annotations to add
#' @param sample_column name of sample name column
#' @param columns_to_add which annotation columns to add, defaults to all
#' @return Sets sample_annotations slot {no return}
#' @examples
#' annotations <- readRDS(system.file("testdata", "dt.rds", package = "BAGEL"))
#' bay <- new("bagel")
#' init_sample_annotations(bay)
#' @export
add_sample_annotations <- function(bay, annotations, sample_column =
                                     "Sample_ID", columns_to_add =
                                     colnames(annotations)) {
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
#' get_sample_annotations(bay)
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
#' get_sample_names(bay)
#' @export
get_sample_names <- function(bay) {
  return(unique(bay@variants$Tumor_Sample_Barcode))
}

#' Creates a new bagel subsetted to only samples with enough variants
#'
#' @param bay Input bagel
#' @param table_name Name of table used for subsetting
#' @param num_counts Minimum sum count value to drop samples
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' subset_bagel_by_counts(bay, "SNV96", 5)
#' @export
subset_bagel_by_counts <- function(bay, table_name, num_counts) {
  tab <- extract_count_table(bay, table_name)
  min_samples <- colnames(tab)[which(colSums(tab) >= num_counts)]

  #Subset all tables
  table_names <- names(bay@count_tables@table_name)
  for (name in table_names) {
    sub_tab <- bay@count_tables@table_list[[name]]
    sub_tab <- sub_tab[, which(colnames(sub_tab) %in% min_samples)]
    bay@count_tables@table_list[[name]] <- sub_tab
  }

  #Subset variants
  bay@variants <- bay@variants[which(bay@variants$Tumor_Sample_Barcode %in%
                                      min_samples), ]

  #Subset sample annotations
  if (nrow(bay@sample_annotations) != 0) {
    bay@sample_annotations <- bay@sample_annotations[which(
      bay@sample_annotations$Samples %in% min_samples), ]
  }
  return(bay)
}

#' Creates a new bagel subsetted to only one value of a sample annotation
#'
#' @param bay Input bagel
#' @param annot_col Annotation class to use for subsetting
#' @param annot_name Annotational value to subset to
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' subset_bagel_by_annotation(bay, "Tumor_Type", "Lung")
#' @export
subset_bagel_by_annotation <- function(bay, annot_col, annot_name) {
  if (!annot_col %in% colnames(bay@sample_annotations)) {
    stop(paste(annot_col, " not found in annotation columns, please review.",
               sep = ""))
  }
  annotation_index <- which(bay@sample_annotations[[which(colnames(
    bay@sample_annotations) %in% annot_col)]] == annot_name)
  if (length(annotation_index) == 0) {
    stop(paste(annot_name, " not present in ", annot_col,
               " column, please review.", sep = ""))
  }
  bay@sample_annotations <- bay@sample_annotations[annotation_index, ]
  annotation_samples <- bay@sample_annotations$"Samples"
  bay@counts_table <- bay@counts_table[, which(colnames(bay@counts_table) %in%
                                              annotation_samples)]
  bay@variants <- bay@variants[which(bay@variants$Tumor_Sample_Barcode %in%
                                      annotation_samples), ]
  return(bay)
}

# Result object/methods -------------------------------

#' Object containing deconvolved/predicted signatures, sample weights, and
#' the bagel object the result was generated from
#'
#' @slot signatures A matrix of signatures by mutational motifs
#' @slot samples A matrix of samples by signature weights
#' @slot type Describes how the signatures/weights were generated
#' @slot bagel The bagel object the results were generated from
#' @slot log_lik Posterior likelihood of the result (LDA only)
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
