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
    missing or malformed, please run create tables e.g. [create_snv96_table]
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
