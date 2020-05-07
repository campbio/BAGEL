#' Uses a genome object to find context and generate standard SNV96 tables
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

#' Uses a genome object to find context and generate standard SNV192 table
#' using transcript strand
#'
#' @param bay Input samples
#' @param g Genome object used for finding variant context
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' g <- select_genome("38")
#' annotate_transcript_strand(bay, "19")
#' create_snv192_table(bay, g)
#' @export
create_snv192_table <- function(bay, g) {
  dat <- bay@variants
  dat <- drop_na_variants(dat, "Transcript_Strand")

  #Fix Chromosomes
  chr <- dat$Chromosome
  tryCatch(
    GenomeInfoDb::seqlevelsStyle(chr) <- "UCSC",
    error = function(e) {
      warning("found no sequence renaming map compatible with seqname",
              " style 'UCSC' for the input reference ", (g@pkgname))
    }
  )

  #Mutation Context
  range_start <- dat$Start_Position
  range_end <- dat$End_Position
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 1, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 1,
                                      as.character = TRUE)
  ref_context <- paste(lflank, dat$Tumor_Seq_Allele1, rflank, sep = "")

  final_mut_type <- rep(NA, nrow(dat))
  final_mut_context <- rep(NA, nrow(dat))

  ## Get mutation type
  initial_maf_type <- paste(dat$Tumor_Seq_Allele1, ">", dat$Tumor_Seq_Allele2,
                           sep = "")

  ## Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- dat$Variant_Type == "SNV" & initial_maf_type %in% forward_change

  final_mut_type[ind] <- initial_maf_type[ind]
  final_mut_context[ind] <- ref_context[ind]

  ## Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- dat$Variant_Type == "SNV" & initial_maf_type %in% rev_change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev_context <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    ref_context[ind]))
  rev_refbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$Tumor_Seq_Allele1[ind]))
  rev_altbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$Tumor_Seq_Allele2[ind]))

  final_mut_type[ind] <- paste(as.character(rev_refbase), ">",
                              as.character(rev_altbase), sep = "")
  final_mut_context[ind] <- rev_context

  maf_mut_id <- paste(final_mut_type, final_mut_context, dat$Transcript_Strand,
                     sep = "_")
  tumor_id <- as.factor(dat$Tumor_Sample_Barcode)

  ## Define all mutation types for 196 substitution scheme
  b1 <- rep(rep(c("A", "C", "G", "T"), each = 24), 2)
  b2 <-  rep(rep(c("C", "T"), each = 12), 8)
  b3 <- rep(c("A", "C", "G", "T"), 48)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut_type <- rep(rep(rep(forward_change, each = 4), 4), 2)
  mut_strand <- rep(c("T", "U"), each = 96)

  mut_id <- apply(cbind(mut_type, mut_trinuc, mut_strand), 1, paste,
                 collapse = "_")

  mutation <- factor(maf_mut_id, levels = mut_id)

  mut_table <- xtabs(~ mutation + tumor_id)

  #Convert to table by dropping xtabs class and call
  attr(mut_table, "call") <- NULL
  attr(mut_table, "class") <- NULL
  tab <- create_count_table(
    bay = bay, table = mut_table, name = "SNV192", description =
      paste("Single Nucleotide Variant table with one base upstream and",
            " downstream and transcript strand", sep = ""),
    return_instead = TRUE)
  eval.parent(substitute(bay@count_tables <- tab))
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

  ref <- dbs$Tumor_Seq_Allele1
  alt <- dbs$Tumor_Seq_Allele2

  #Reverse Complement broad categories to the other strand
  rc_ref <- which(ref %in% c("GT", "GG", "AG", "GA", "CA", "AA"))
  ref[rc_ref] <- rc(ref[rc_ref])
  alt[rc_ref] <- rc(alt[rc_ref])

  #Reverse complement more specific categories to the other strand
  rc_alt <- NULL
  rc_ref_AT <- which(ref == "AT")
  rc_alt <- c(rc_alt, rc_ref_AT[which(alt[rc_ref_AT] %in% c("TG", "GG", "TC"))])

  rc_ref_CG <- which(ref == "CG")
  rc_alt <- c(rc_alt, rc_ref_CG[which(alt[rc_ref_CG] %in% c("AC", "GA", "AA"))])

  rc_ref_GC <- which(ref == "GC")
  rc_alt <- c(rc_alt, rc_ref_GC[which(alt[rc_ref_GC] %in% c("TT", "CT", "TG"))])

  rc_ref_TA <- which(ref == "TA")
  rc_alt <- c(rc_alt, rc_ref_TA[which(alt[rc_ref_TA] %in% c("AG", "CC", "AC"))])

  if (length(rc_alt) > 0) {
    alt[rc_alt] <- rc(alt[rc_alt])
  }

  full <- paste(ref, ">NN_", alt, sep = "")
  full_motif <- c(paste0("AC>NN", "_", c("CA", "CG", "CT", "GA", "GG", "GT",
                                         "TA", "TG", "TT")),
            paste0("AT>NN", "_", c("CA", "CC", "CG", "GA", "GC", "TA")),
            paste0("CC>NN", "_", c("AA", "AG", "AT", "GA", "GG", "GT", "TA",
                                   "TG", "TT")),
            paste0("CG>NN", "_", c("AT", "GC", "GT", "TA", "TC", "TT")),
            paste0("CT>NN", "_", c("AA", "AC", "AG", "GA", "GC", "GG", "TA",
                                   "TC", "TG")),
            paste0("GC>NN", "_", c("AA", "AG", "AT", "CA", "CG", "TA")),
            paste0("TA>NN", "_", c("AT", "CG", "CT", "GC", "GG", "GT")),
            paste0("TC>NN", "_", c("AA", "AG", "AT", "CA", "CG", "CT", "GA",
                                   "GG", "GT")),
            paste0("TG>NN", "_", c("AA", "AC", "AT", "CA", "CC", "CT", "GA",
                                   "GC", "GT")),
            paste0("TT>NN", "_", c("AA", "AC", "AG", "CA", "CC", "CG", "GA",
                                   "GC", "GG")))

  sample_names <- unique(dbs$Tumor_Sample_Barcode)
  num_samples <- length(sample_names)
  variant_tables <- vector("list", length = num_samples)
  for (i in seq_len(num_samples)) {
    sample_index <- which(dbs$Tumor_Sample_Barcode == sample_names[i])
    variant_tables[[i]] <- table(factor(full[sample_index],
                                        levels = full_motif))
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

#' Reverse complement of a string using biostrings
#'
#' @param dna Input DNA string
#' @examples
#' rc("ATGC")
#' @export
rc <- function(dna) {
  if (class(dna) == "character" && length(dna) == 1) {
    rev_com <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(dna)))
  } else if (class(dna) == "character" && length(dna) > 1) {
    rev_com <- sapply(dna, rc)
    names(rev_com) <- NULL
  } else {
    stop("Must be character or character vector")
  }
  return(rev_com)
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
