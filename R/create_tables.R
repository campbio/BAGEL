#' Uses a genome object to find context and add it to the variant table
#'
#' @param bay Input samples
#' @param g Genome object used for finding variant context
#' @param flank_start Start of flank area to add, can be positive or negative
#' @param flank_end End of flank area to add, can be positive or negative
#' @param build_table Automatically build a table using the annotation and add
#' it to the bagel
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' g <- select_genome("38")
#' add_flank_to_variants(bay, g, 1, 2)
#' add_flank_to_variants(bay, g, -2, -1)
#' @export
add_flank_to_variants <- function(bay, g, flank_start, flank_end,
                                  build_table = TRUE) {
  stopifnot(sign(flank_start) == sign(flank_end), flank_start < flank_end)

  direction <- ifelse(sign(flank_start) == 1, "r", "l")

  #Determine output and calculations based on selected context area
  output_column <- paste(direction, "flank_", abs(flank_start), "_to_",
                         abs(flank_end), sep = "")

  dat <- bay@variants
  mut_type <- paste(dat$Tumor_Seq_Allele1, ">", dat$Tumor_Seq_Allele2, sep = "")
  chr <- dat$Chromosome

  if (sign(flank_start) == 1) {
    center <- dat$End_Position
  } else {
    center <- dat$Start_Position
  }
  ref <- dat$Tumor_Seq_Allele1
  type <- mut_type

  #Mutation Context
  flank <- VariantAnnotation::getSeq(g, chr, center + flank_start,
                                     center + flank_end, as.character = TRUE)
  final_mut_context <- rep(NA, length(ref))

  # Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- type %in% forward_change

  final_mut_context[ind] <- flank[ind]

  # Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- type %in% rev_change

  # Reverse complement the context so only 6 mutation categories instead of 12
  rev_flank <- flank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()

  final_mut_context[ind] <- as.character(rev_flank)
  dat[[output_column]] <- final_mut_context
  eval.parent(substitute(bay@variants <- dat))
  if (build_table) {
    dat_bagel = methods::new("bagel", variants = dat, count_tables =
                               bay@count_tables,
                             sample_annotations = bay@sample_annotations)
    tab <- create_variant_table(dat_bagel, variant_annotation = output_column,
                         name = output_column, return_instead = FALSE)
    eval.parent(substitute(bay@count_tables <- tab))
  }
}

#' Adds an annotation to the variant table with length of each variant
#'
#' @param bay Input samples
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_variant_length(bay)
#' @export
annotate_variant_length <- function(bay) {
  dat <- bay@variants
  var_length <- rep(NA, nrow(dat))
  var_length[which(dat$Variant_Type == "SNV")] <- 1
  var_length[which(dat$Variant_Type == "DBS")] <- 2
  indels <- which(dat$Variant_Type == "indel")
  var_length[indels] <- nchar(dat$Tumor_Seq_Allele2[indels]) -
    nchar(dat$Tumor_Seq_Allele1[indels])
  dat[["Variant_Length"]] <- var_length
  eval.parent(substitute(bay@variants <- dat))
}

#' Drops a column from the variant table that the user no longer needs
#'
#' @param bay Input samples
#' @param column_name Name of column to drop
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' drop_annotation(bay, "Chromosome")
#' @export
drop_annotation <- function(bay, column_name) {
  dat <- bay@variants
  stopifnot(column_name %in% colnames(dat))
  data.table::set(dat, j = column_name, value = NULL)
  eval.parent(substitute(bay@variants <- dat))
}

#' Uses a genome object to find context and generate tables for input samples
#'
#' @param bay Input samples
#' @param g Genome object used for finding variant context
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' g <- select_genome("38")
#' create_snv96_table(bay, g)
#' @export
create_snv96_table <- function(bay, g) {
  dat <- bay@variants

  mut_type <- paste(dat$Tumor_Seq_Allele1, ">", dat$Tumor_Seq_Allele2, sep = "")

  #Fix Chromosomes
  #### TODO

  chr <- dat$Chromosome
  tryCatch(
    GenomeInfoDb::seqlevelsStyle(chr) <- "UCSC",
    error = function(e) {
      warning("found no sequence renaming map compatible with seqname",
              " style 'UCSC' for the input reference ", (g@pkgname))
    }
  )

  range_start <- dat$Start_Position
  range_end <- dat$End_Position
  ref <- dat$Tumor_Seq_Allele1
  alt <- dat$Tumor_Seq_Allele2
  type <- mut_type

  #Mutation Context
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 1, range_start - 1,
                   as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 1,
                   as.character = TRUE)

  final_mut_type <- rep(NA, length(ref))
  final_mut_context <- rep(NA, length(ref))

  # Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- type %in% forward_change

  final_mut_type[ind] <- paste(as.character(ref[ind]), ">",
                               as.character(alt[ind]), sep = "")
  final_mut_context[ind] <- paste(lflank[ind], ref[ind], rflank[ind], sep = "")

  # Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- type %in% rev_change

  # Reverse complement the context so only 6 mutation categories instead of 12
  rev_refbase <- ref[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()
  rev_altbase <- alt[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()
  rev_lflank <- lflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()
  rev_rflank <- rflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()

  final_mut_lflank <- as.character(rev_lflank)
  final_mut_rflank <- as.character(rev_rflank)

  final_mut_type[ind] <- paste(as.character(rev_refbase), ">",
                               as.character(rev_altbase), sep = "")
  final_mut_context[ind] <- paste(final_mut_lflank, rev_refbase,
                                  final_mut_rflank, sep = "")
  maf_mut_id <- paste(final_mut_type, final_mut_context, sep = "_")


  ###### Now we separate into samples
  ## Define all mutation types for 96 substitution scheme
  b1 <- rep(c("A", "C", "G", "T"), each = 24)
  b2 <- rep(rep(c("C", "T"), each = 12), 4)
  b3 <- rep(c("A", "C", "G", "T"), 24)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut_type <- rep(rep(forward_change, each = 4), 4)

  sample_names <- unique(dat$Tumor_Sample_Barcode)
  num_samples <- length(sample_names)
  maf_mut_summaries <- vector("list", length = num_samples)
  for (i in seq_len(num_samples)) {
    sample_index <- which(dat$Tumor_Sample_Barcode == sample_names[i])
    mut_id <- apply(cbind(mut_type, mut_trinuc), 1, paste, collapse = "_")
    mutation <- factor(maf_mut_id[sample_index], levels = mut_id)
    mut_summary <- data.frame(mutation, Type = final_mut_type[sample_index],
                              Context = final_mut_context[sample_index],
                              stringsAsFactors = FALSE)
    maf_mut_summaries[[i]] <- mut_summary
  }
  mut_summary_mat <- do.call(cbind, lapply(maf_mut_summaries, function(x)
    table(x[, "mutation"])))
  colnames(mut_summary_mat) <- sample_names
  tab <- create_count_table(bay = bay, table = mut_summary_mat, name = "SNV96",
                     description = paste("Single Nucleotide Variant table with",
                     " one base upstream and downstream",
                                         sep = ""), return_instead = TRUE)
  eval.parent(substitute(bay@count_tables <- tab))
}

#' Generates a variant type table
#'
#' @param tab Input variant table
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' BAGEL:::add_variant_type(bay@variants)
#' @keywords internal
add_variant_type <- function(tab) {
  type <- rep(NA, nrow(tab))
  type[which(nchar(tab$Tumor_Seq_Allele1) == 1 &
               nchar(tab$Tumor_Seq_Allele2) == 1)] <- "SNV"
  type[which(nchar(tab$Tumor_Seq_Allele1) == 2 &
               nchar(tab$Tumor_Seq_Allele2) == 2)] <- "DBS"
  type[which(tab$Tumor_Seq_Allele1 == "-" |
               tab$Tumor_Seq_Allele2 == "-")] <- "indel"
  type[which(is.na(type))] <- "indel"
  tab$Variant_Type <- type
  return(tab)
}

#' Annotate variants with variant type ("SNV", "INS", "DEl", "DBS")
#'
#' @param bay Input bagel
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_variant_type(bay)
#' @export
annotate_variant_type <- function(bay) {
  type_added <- add_variant_type(bay@variants)
  eval.parent(substitute(bay@variants <- type_added))
}

#' Subsets a variant table based on Variant Type
#'
#' @param tab Input variant table
#' @param type Variant type to return e.g. "SNV", "INS", "DEL", "DBS"
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_variant_type(bay)
#' subset_variant_by_type(get_variants(bay), "SNV")
#' @export
subset_variant_by_type <- function(tab, type) {
  if(!"Variant_Type" %in% colnames(tab)) {
    stop(paste("No Variant_Type annotation found, ",
               "please run annotate_variant_type first."))
  }
  if(!any(tab$Variant_Type %in% type)) {
    stop(paste("No variants of type: ", type))
  }
  return(tab[which(tab$Variant_Type == type), ])
}

#' Creates and adds a table for standard doublet base subsitution (DBS)
#'
#' @param bay Input bagel
#' @examples
#' bay <- readRDS(system.file("testdata", "dbs_bagel.rds", package = "BAGEL"))
#' create_dbs_table(bay)
#' @export
create_dbs_table <- function(bay) {
  dbs <- subset_variant_by_type(bay@variants, "DBS")

  dbs_motifs <- cbind(paste(dbs$Tumor_Seq_Allele1, dbs$Tumor_Seq_Allele2, sep=">"))
  default_factor <- levels(factor(dbs_motifs))
  sample_names <- unique(dbs$Tumor_Sample_Barcode)
  num_samples <- length(sample_names)
  variant_tables <- vector("list", length = num_samples)
  for (i in seq_len(num_samples)) {
    sample_index <- which(dbs$Tumor_Sample_Barcode == sample_names[i])
    variant_tables[[i]] <- table(factor(dbs_motifs[sample_index],
                                        levels = default_factor))
  }
  table <- do.call(cbind, variant_tables)
  colnames(table) <- sample_names

  tab <- create_count_table(bay = bay, table = table, name = "DBS",
                            description = paste("Standard count table for ",
                                                "double-base substitutions",
                                                sep = ""),
                            return_instead = TRUE)
  eval.parent(substitute(bay@count_tables <- tab))
}

create_indel_table <- function(bay, g) {
  temp <- methods::new("bagel", variants = subset_variant_by_type(bay@variants,
                                                                  "indel"),
                       count_tables = bay@count_tables)
  tab <- create_variant_table(bay = temp, variant_annotation = "Variant_Length",
                              name = "indel", description =
                                "Standard count table for indels",
                              return_instead = TRUE)
  eval.parent(substitute(bay@count_tables <- tab))
}

#' Add strand annotation to SNV variants
#'
#' @param bay Input bagel
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_snv_strand(bay)
#' @export
annotate_snv_strand <- function(bay) {
  dat <- bay@variants
  snvs <- which(dat$Variant_Type == "SNV")
  motifs <- paste(dat$Tumor_Seq_Allele1,
                  dat$Tumor_Seq_Allele2, sep = ">")
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  strand <- rep(NA, nrow(dat))
  strand[which(motifs %in% forward_change)] <- "+"
  strand[which(motifs %in% rev_change)] <- "-"
  dat[["SNV_Strand"]] <- strand
  eval.parent(substitute(bay@variants <- dat))
}
