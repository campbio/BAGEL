#' Uses a genome object to find context and generate standard SBS96 tables
#'
#' @param bay Input samples
create_sbs96_table <- function(bay) {
  dat <- subset_variant_by_type(bay@variants, type = "SBS")
  g <- bay@genome
  ref <- as.character(dat$ref)
  alt <- as.character(dat$alt)
  mut_type <- paste(ref, ">", alt, sep = "")

  #Mutation Context
  context <- BSgenome::getSeq(x = g, names = as.character(dat$chr),
                              start = dat$start - 1,
                              end = dat$end + 1,
                              as.character = TRUE)

  final_mut_type <- rep(NA, length(ref))
  final_mut_context <- rep(NA, length(ref))

  # Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- mut_type %in% forward_change
  final_mut_type[ind] <- paste(as.character(ref[ind]), ">",
                               as.character(alt[ind]), sep = "")
  final_mut_context[ind] <- context[ind]

  # Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- mut_type %in% rev_change

  # Reverse complement the context so only 6 mutation categories instead of 12
  rev_context <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(context[ind])))
  rev_refbase <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(ref[ind])))
  rev_altbase <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(alt[ind])))

  final_mut_type[ind] <- paste0(rev_refbase, ">", rev_altbase)
  final_mut_context[ind] <- rev_context
  final_motif <- paste0(final_mut_type, "_", final_mut_context)

  ###### Now we separate into samples
  ## Define all mutation types for 96 substitution scheme
  b1 <- rep(rep(c("A", "C", "G", "T"), each=4), 6)
  b2 <- rep(c("C", "T"), each = 48)
  b3 <- rep(c("A", "C", "G", "T"), 24)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut <- rep(forward_change, each = 16)
  annotation <- data.frame("motif" = paste0(mut, "_", mut_trinuc),
                           "mutation" = mut,
                           "context" = mut_trinuc)
  rownames(annotation) <- annotation$motif
  

  mut_id <- apply(cbind(mut, mut_trinuc), 1, paste,
                  collapse = "_")
  mutation <- factor(final_motif, levels = annotation$motif)

  mut_table <- as.matrix(as.data.frame.matrix(xtabs(~ mutation + dat$sample)))
  #xtabs adds dat$sample as dimname[2] so we remove it
  dimnames(mut_table) <- list(rownames(mut_table), colnames(mut_table))

  # Use COSMIC color scheme
  color_mapping <- c("C>A" = "#5ABCEBFF",
                     "C>G" = "#050708FF",
                     "C>T" = "#D33C32FF",
                     "T>A" = "#CBCACBFF",
                     "T>C" = "#ABCD72FF",
                     "T>G" = "#E7C9C6FF")
 
# Need to think about this more carefully - what happens if indel are zero but
#  not SNV and the user wants to combine them? 
#  zero_samps <- which(colSums(mut_table) == 0)
#  if (length(zero_samps) > 0) {
#    warning(paste0("Dropping the following zero count samples: ",
#                   paste(names(zero_samps), collapse = ", ")))
#    mut_table <- mut_table[, -zero_samps, drop = FALSE]
#  }
  tab <- .create_count_table(bay = bay, name = "SBS96",
                            count_table = mut_table, 
                            annotation = annotation,
                            features = data.frame(mutation = final_motif),
                            type = as.character(dat$Variant_Type),
                            color_variable = "mutation",
                            color_mapping = color_mapping,
                     description = paste0("Single Base Substitution table with",
                     " one base upstream and downstream"))
  return(tab)
}

#' Uses a genome object to find context and generate standard SBS192 table
#' using transcript strand
#'
#' @param bay Input samples
#' @param strand_type Transcript_Strand or Replication_Strand
create_sbs192_table <- function(bay, strand_type) {
  if (!strand_type %in% c("Transcript_Strand", "Replication_Strand")) {
    stop("Please select either Transcript_Strand or Replication_Strand")
  }
  g <- bay@genome
  dat <- bay@variants
  dat <- drop_na_variants(dat, strand_type)

  chr <- dat$chr
  range_start <- dat$start
  range_end <- dat$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 1, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 1,
                                      as.character = TRUE)
  ref_context <- paste(lflank, dat$ref, rflank, sep = "")

  final_mut_type <- rep(NA, nrow(dat))
  final_mut_context <- rep(NA, nrow(dat))

  ## Get mutation type
  initial_maf_type <- paste(dat$ref, ">", dat$alt,
                           sep = "")

  ## Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- dat$Variant_Type == "SBS" & initial_maf_type %in% forward_change

  final_mut_type[ind] <- initial_maf_type[ind]
  final_mut_context[ind] <- ref_context[ind]

  ## Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- dat$Variant_Type == "SBS" & initial_maf_type %in% rev_change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev_context <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    ref_context[ind]))
  rev_refbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$ref[ind]))
  rev_altbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$alt[ind]))

  final_mut_type[ind] <- paste(as.character(rev_refbase), ">",
                              as.character(rev_altbase), sep = "")
  final_mut_context[ind] <- rev_context

  maf_mut_id <- paste(final_mut_type, final_mut_context, dat[[strand_type]],
                     sep = "_")
  tumor_id <- as.factor(dat$sample)

  ## Define all mutation types for 196 substitution scheme
  b1 <- rep(rep(c("A", "C", "G", "T"), each = 24), 2)
  b2 <-  rep(rep(c("C", "T"), each = 12), 8)
  b3 <- rep(c("A", "C", "G", "T"), 48)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut_type <- rep(rep(rep(forward_change, each = 4), 4), 2)
  if (strand_type == "Transcript_Strand") {
    mut_strand <- rep(c("T", "U"), each = 96)
  } else if (strand_type == "Replication_Strand") {
    mut_strand <- rep(c("leading", "lagging"), each = 96)
  }
  mut_id <- apply(cbind(mut_type, mut_trinuc, mut_strand), 1, paste,
                 collapse = "_")

  mutation <- factor(maf_mut_id, levels = mut_id)

  mut_table <- xtabs(~ mutation + tumor_id)

  #Convert to table by dropping xtabs class and call
  attr(mut_table, "call") <- NULL
  attr(mut_table, "class") <- NULL
  zero_samps <- which(colSums(mut_table) == 0)
  if (length(zero_samps) > 0) {
    warning(paste0("Dropping the following zero count samples: ",
                   paste(colnames(mut_table[zero_samps]), sep = ", ")))
    mut_table <- mut_table[, -zero_samps]
  }
  tab <- create_count_table(
    bay = bay, table = mut_table, name = paste0("SBS192_", ifelse(
      strand_type == "Transcript_Strand", "Trans", "Rep")), description =
      paste("Single Base Substitution table with one base upstream and",
            " downstream and transcript strand", sep = ""),
    return_instead = TRUE)
  return(tab)
}

#' Creates and adds a table for standard doublet base subsitution (DBS)
#'
#' @param bay Input bagel
create_dbs_table <- function(bay) {
  dbs <- subset_variant_by_type(bay@variants, "DBS")

  ref <- dbs$ref
  alt <- dbs$alt

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

  sample_names <- unique(dbs$sample)
  num_samples <- length(sample_names)
  variant_tables <- vector("list", length = num_samples)
  for (i in seq_len(num_samples)) {
    sample_index <- which(dbs$sample == sample_names[i])
    variant_tables[[i]] <- table(factor(full[sample_index],
                                        levels = full_motif))
  }
  table <- do.call(cbind, variant_tables)
  colnames(table) <- sample_names
  zero_samps <- which(colSums(table) == 0)
  if (length(zero_samps) > 0) {
    warning(paste0("Dropping the following zero count samples: ",
                   paste(colnames(table[zero_samps]), sep = ", ")))
    table <- table[, -zero_samps]
  }
  tab <- create_count_table(bay = bay, table = table, name = "DBS",
                            description = paste("Standard count table for ",
                                                "double base substitutions",
                                                sep = ""),
                            return_instead = TRUE)
  return(tab)
}

#' Reverse complement of a string using biostrings
#'
#' @param dna Input DNA string
#' @examples
#' rc("ATGC")
#' @export
rc <- function(dna) {
  if (class(dna) == "character" && length(dna) == 1) {
    rev_com <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAString(dna)))
  } else if (class(dna) == "character" && length(dna) > 1) {
    rev_com <- sapply(dna, rc)
    names(rev_com) <- NULL
  } else {
    stop("Must be character or character vector")
  }
  return(rev_com)
}

#' Builds a standard table from user variants
#'
#' @param bay Input samples
#' @param table_name Name of standard table to build SBS96, SBS192, DBS
#' @param strand_type Only for SBS192 Transcript_Strand or Replication_Strand
#' Indel
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' build_standard_table(bay, "SBS96")
#'
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_transcript_strand(bay, "19")
#' build_standard_table(bay, "SBS192", "Transcript_Strand")
#'
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_replication_strand(bay, BAGEL::rep_range)
#' build_standard_table(bay, "SBS192", "Replication_Strand")
#'
#' bay <- readRDS(system.file("testdata", "dbs_bagel.rds",
#' package = "BAGEL"))
#' build_standard_table(bay, table_name = "DBS")
#'
#' bay <- readRDS(system.file("testdata", "indel_bagel.rds", package = "BAGEL"))
#' build_standard_table(bay, table_name = "INDEL")
#' @export
build_standard_table <- function(bay, table_name, strand_type = NA) {
  if (table_name %in% c("SNV96", "SNV", "96", "SBS", "SBS96")) {
    tab <- create_sbs96_table(bay)
  } else if (table_name %in% c("SBS192", "192")) {
    tab <- create_sbs192_table(bay, strand_type)
  } else if (table_name %in% c("DBS", "doublet")) {
    tab <- create_dbs_table(bay)
  } else if (table_name %in% c("INDEL", "IND", "indel", "Indel")) {
    tab <- create_indel_table(bay)
  } else {
    stop(paste0("There is no standard table named: ", table_name,
               " please select from SBS96, SBS192, DBS, Indel."))
  }
  eval.parent(substitute(bay@count_tables <- tab))
}

create_indel_table <- function(bay) {
  var <- bay@variants
  g <- bay@genome
  all_ins <- subset_variant_by_type(var, "INS")
  all_del <- subset_variant_by_type(var, "DEL")
  samples <- unique(c(as.character(all_ins$sample),
                      as.character(all_del$sample)))
  dimlist <- list(row_names = c(.get_indel_motifs("bp1", 0, 0),
                                .get_indel_motifs("bp1", 1, 0),
                                .get_indel_motifs("del", NA, NA),
                                .get_indel_motifs("ins", NA, NA),
                                .get_indel_motifs("micro", NA, NA)),
                  column_names = samples)
  table <- matrix(NA, nrow = 83, ncol = length(samples), dimnames = dimlist)
  for(sample in samples) {
    ins <- all_ins[which(all_ins$sample == sample), ]
    del <- all_del[which(all_del$sample == sample), ]

    ins_len = nchar(ins$alt)
    del_len = nchar(del$ref)
    ins1 = ins[which(ins_len == 1), ]
    ins2 = ins[which(ins_len > 1), ]

    del1 = del[which(del_len == 1), ]
    del2 = del[which(del_len > 1), ]

    if(nrow(del1) == 0) {
      del1_counts <- setNames(rep(0, 12), .get_indel_motifs("bp1", 0, 0))
    } else {
      del1_counts <- .count1(del1, del1$ref, ins = FALSE, g = g)
    }

    if(nrow(ins1) == 0) {
      ins1_counts <- setNames(rep(0, 12), .get_indel_motifs("bp1", 1, 0))
    } else {
      ins1_counts <- .count1(mut = ins1, type = ins1$alt, ins = TRUE, g = g)
    }

    if(nrow(ins2) == 0) {
      ins2_counts <- setNames(rep(0, 24), .get_indel_motifs("ins", NA, NA))
    } else {
      ins2_counts <- .count2_ins(mut = ins2, type = ins2$alt, g = g)
    }

    if(nrow(del2) == 0) {
      del2_counts <- list(del = setNames(rep(0, 24),
                                         .get_indel_motifs("del", NA, NA)),
                          micro = setNames(rep(0, 11),
                                           .get_indel_motifs("micro", NA, NA)))
    } else {
      del2_counts = .count2_del(mut = del2, type = del2$ref, g)
    }
    table[, sample] <- c(del1_counts, ins1_counts, del2_counts$del, ins2_counts,
                       del2_counts$micro)
  }
  tab <- create_count_table(bay = bay, table = table, name = "INDEL",
                            description = paste("Standard count table for ",
                                                "small insertions and deletions"
                                                ,sep = ""),
                            return_instead = TRUE)
  return(tab)
}

.get_indel_motifs <- function(indel, ins, plus) {
  if(indel == "bp1") {
    return(paste(ifelse(ins, "INS", "DEL"), c("C", "C", "C", "C", "C", "C", "T",
                                              "T", "T", "T", "T", "T"), 1,
                 c(c(0, 1, 2, 3, 4, paste0(5 + plus, "+"), 0, 1, 2, 3, 4,
                     paste0(5 + plus, "+"))), sep = "_"))
  } else if(indel == "ins") {
    return(paste(paste("INS_repeats", c("2", "2", "2", "2", "2", "2", "3", "3", "3", "3",
                                        "3", "3", "4", "4", "4", "4", "4", "4", "5+",
                                        "5+", "5+", "5+", "5+", "5+"),
                       c("0", "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4",
                         "5+", "0", "1", "2", "3", "4", "5+", "0", "1", "2", "3",
                         "4", "5+"), sep = "_")))
  } else if(indel == "micro") {
    return(paste("DEL_MH",  c("2", "3", "3", "4", "4", "4", "5+", "5+",
                              "5+", "5+", "5+"),
                 c("1", "1", "2", "1", "2", "3", "1", "2",
                   "3", "4", "5+"), sep = "_"))
  } else if(indel == "del") {
    return(paste("DEL_repeats", c("2", "2", "2", "2", "2", "2", "3", "3", "3", "3",
                                  "3", "3", "4", "4", "4", "4", "4", "4", "5+",
                                  "5+", "5+", "5+", "5+", "5+"),
                 c("0", "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4",
                   "5+", "0", "1", "2", "3", "4", "5+", "0", "1", "2", "3",
                   "4", "5+"), sep = "_"))
  } else {
    stop("Unrecognized indel type")
  }
}

.count1 <- function(mut, type, ins, g) {
  ifelse(ins, plus <- 0, plus <- 0)
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 10, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 10,
                                      as.character = TRUE)
  ind <- which(type %in% c("A", "G"))
  lflank[ind] <- lflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  rflank[ind] <- rflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  final_type <- type
  final_type[ind] <- final_type[ind] %>% Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>% as.character()
  repeats = rep(NA, length(final_type))
  for(i in 1:length(type)) {
    repeats[i] <- .count_repeat(final_type[i], rflank[i]) +
      .count_repeat(final_type[i], rev(lflank[i])) + plus
  }
  repeats[repeats >= 5 + plus] <- paste0(5 + plus, "+")
  bp1_motif <- .get_indel_motifs("bp1", ins, plus)
  return(table(factor(paste(ifelse(ins, "INS", "DEL"), final_type, 1,
                            repeats, sep = "_"), levels = bp1_motif)))
}

.count2_ins <- function(mut, type, g) {
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 10, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 10,
                                      as.character = TRUE)
  repeats = rep(NA, length(type))
  for(i in 1:length(type)) {
    repeats[i] <- .count_repeat(type[i], rflank[i]) +
      .count_repeat(type[i], rev(lflank[i]))
  }
  repeats[repeats >= 5] <- paste0(5, "+")
  len <- nchar(type)
  len[which(len >= 5)] <- "5+"
  ins_motif <- .get_indel_motifs("ins", NA, NA)
  return(table(factor(paste("INS_repeats", len, repeats, sep = "_"),
                      levels = ins_motif)))
}

.count2_del <- function(mut, type, g) {
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end

  len <- nchar(type)
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - len,
                                      range_start - 1, as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1,
                                      range_end + len, as.character = TRUE)
  has_repeat <- type == lflank | type == rflank
  maybe_micro <- which(!has_repeat)

  micro <- rep(NA, length(type))
  for(i in 1:length(type)) {
    micro[i] <- max(.micro_left(type[i], lflank[i]),
                    .micro_right(type[i], rflank[i]))
  }

  repeats = rep(NA, length(type))
  for(i in 1:length(type)) {
    repeats[i] <- .count_repeat(type[i], rflank[i]) +
      .count_repeat(type[i], rev(lflank[i])) + 1
  }
  micro_ind <- which(repeats == 1 & micro > 0)
  repeat_ind <- which(micro == 0)
  final_micro <- micro[micro_ind]
  final_repeats <- repeats[repeat_ind]
  final_repeats[final_repeats >= 6] <- paste0(6, "+")
  final_len <- len
  final_len[which(final_len >= 5)] <- "5+"
  final_micro[which(final_micro >= 5)] <- "5+"
  micro_motif <- .get_indel_motifs("micro", NA, NA)
  del_motif <- .get_indel_motifs("del", NA, NA)
  del_tab <- table(factor(paste("DEL_repeats", final_len[repeat_ind],
                                final_repeats, sep = "_"), levels = del_motif))
  micro_tab <- table(factor(paste("DEL_MH", final_len[micro_ind], final_micro,
                                  sep = "_"), levels = micro_motif))
  return(list(del = del_tab, micro = micro_tab))
}

.micro_right <-  function(letter, string) {
  len <- nchar(letter)
  elem <- len - 1
  while(elem > 0) {
    if(substr(letter, 1, elem) == substr(string, 1, elem)) {
      return(elem)
    } else {
      elem <- elem - 1
    }
  }
  return(0)
}

.micro_left <-  function(letter, string) {
  len <- nchar(letter)
  elem <- 1
  while(elem < len) {
    if(substr(letter, elem, len) == substr(string, elem, len)) {
      return(elem)
    } else {
      elem <- elem + 1
    }
  }
  return(0)
}

.count_repeat <- function(letter, string) {
  len <- nchar(letter)
  if(letter != substr(string, 1, len)) {
    return(0)
  } else {
    next_matches = TRUE
    counts <- 1
    while(next_matches) {
      #print(paste0("Counts: ", counts))
      if(letter != substr(string, counts*len + 1, counts*len + len)) {
        return(counts)
      } else {
        counts <- counts + 1
      }
    }
  }
}
