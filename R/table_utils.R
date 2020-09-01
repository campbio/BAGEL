#' Generates a 96 motif table based on input counts for plotting
#'
#' @param sample_df Input counts table
#' @return Returns a 96 motif summary table
table_96 <- function(sample_df) {
  motif <- names(sample_df)
  expanded <- rep(motif, sample_df)
  context <- substr(expanded, 5, 7)
  final_mut_type <- substr(expanded, 1, 3)
  final_mut_context <- context

  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  ## Define all mutation types for 96 substitution scheme
  b1 <- rep(c("A", "C", "G", "T"), each = 24)
  b2 <- rep(rep(c("C", "T"), each = 12), 4)
  b3 <- rep(c("A", "C", "G", "T"), 24)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut_type <- rep(rep(forward_change, each = 4), 4)

  mut_id <- apply(cbind(mut_type, mut_trinuc), 1, paste, collapse = "_")
  expanded <- rep(motif, sample_df)
  mutation <- factor(expanded, levels = mut_id)

  mut_summary <- data.frame(mutation, Type = final_mut_type,
                            Context = final_mut_context,
                            stringsAsFactors = FALSE)
  return(mut_summary)
}

extract_count_table <- function(bagel, table_name) {
  #Check that object is a bagel
  if (!methods::is(bagel, "bagel")) {
    stop(strwrap(prefix = " ", initial = "", "The input object is not a
    bagel object, please use new('bagel') to create one."))
  }

  tabs <- bagel@count_tables

  #Check that at least one table exists
  if (length(tabs@table_name) == 0) {
    stop(strwrap(prefix = " ", initial = "", "The counts table is either
    missing or malformed, please run create tables e.g. [build_standard_table]
    prior to this function."))
  }

  #Check that table exists within this bagel
  if (!table_name %in% tabs@table_name) {
    stop(paste(table_name, " does not exist. Current table names are: ",
               tabs@table_name, sep = ""))
  }

  counts_table <- tabs@table_list[[table_name]]
  return(counts_table)
}

subset_count_tables <- function(bay, samples) {
  tables <- bay@count_tables
  table_names <- names(tables@table_name)
  for (name in table_names) {
    sub_tab <- tables@table_list[[name]]
    sub_tab <- sub_tab[, which(colnames(sub_tab) %in% samples)]
    tables@table_list[[name]] <- sub_tab
  }
  return(tables)
}

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
    stop("Error: Table must be a matrix")
  }

  tab@table_name[[name]] <- name
  tab@description[[name]] <- description
  if (return_instead) {
    return(tab)
  } else {
    eval.parent(substitute(bay@count_tables <- tab))
  }
}

#' Builds a custom table from specified user variants
#'
#' @param bay Input samples
#' @param variant_annotation User column to use for building table
#' @param name Table name to refer to (must be unique)
#' @param description Optional description of the table content
#' @param data_factor Full set of table values, in case some are missing from
#' the data. If NA, a superset of all available unique data values will be used
#' @param return_instead Instead of adding to bagel object, return the created
#' table
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' g <- select_genome("38")
#' annotate_transcript_strand(bay, "19", build_table = FALSE)
#' build_custom_table(bay, "Transcript_Strand", "Transcript_Strand",
#' data_factor = factor(c("T", "U")))
#' @export
build_custom_table <- function(bay, variant_annotation, name,
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
    sample_names <- unique(variants$sample)
    num_samples <- length(sample_names)
    default_factor <- levels(factor(column_data))
    variant_tables <- vector("list", length = num_samples)
    for (i in seq_len(num_samples)) {
      sample_index <- which(variants$sample == sample_names[i])
      if (!all(is.na(data_factor))) {
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
    for (i in seq_len(to_comb)) {
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
