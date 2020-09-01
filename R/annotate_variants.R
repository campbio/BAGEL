#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom S4Vectors queryHits decode subjectHits
#' @importFrom GenomicRanges strand
NULL

#' Uses a genome object to find context and add it to the variant table
#'
#' @param bay Input samples
#' @param g Genome object used for finding variant context
#' @param flank_start Start of flank area to add, can be positive or negative
#' @param flank_end End of flank area to add, can be positive or negative
#' @param build_table Automatically build a table using the annotation and add
#' it to the bagel
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel_snv96_tiny.rds",
#' package = "BAGEL"))
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
    dat_bagel <- methods::new("bagel", variants = dat, count_tables =
                               bay@count_tables,
                             sample_annotations = bay@sample_annotations)
    tab <- build_custom_table(dat_bagel, variant_annotation = output_column,
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
  type[which(tab$Tumor_Seq_Allele1 == "-")] <- "INS"
  type[which(tab$Tumor_Seq_Allele2 == "-")] <- "DEL"
  type[which(is.na(type))] <- "unknown"
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
  if (!"Variant_Type" %in% colnames(tab)) {
    stop(paste("No Variant_Type annotation found, ",
               "please run annotate_variant_type first."))
  }
  if (!any(tab$Variant_Type %in% type)) {
    stop(paste("No variants of type: ", type))
  }
  return(tab[which(tab$Variant_Type == type), ])
}

#' Add transcript strand annotation to SNV variants (defined in genes only)
#'
#' @param bay Input bagel
#' @param genome_build Which genome build to use: hg19, hg38, or a custom TxDb
#' object
#' @param build_table Automatically build a table from this annotation
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_transcript_strand(bay, 19)
#' @export
annotate_transcript_strand <- function(bay, genome_build, build_table = TRUE) {
  if (genome_build %in% c("19", "hg19")) {
    genes <- genes(
      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  } else if (genome_build %in% c("38", "hg38")) {
    genes <- genes(
      TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
  } else if (methods::isClass(genome_build, "TxDb")) {
    genes <- genome_build
  } else {
    stop("Please select hg19, hg38, or provde a TxDb object")
  }

  dat <- bay@variants
  snv_index <- which(dat$Variant_Type == "SNV")
  snvs <- subset_variant_by_type(dat, "SNV")

  #Create VRanges object to determine strand of variants within genes
  vrange <- VariantAnnotation::VRanges(seqnames = snvs$Chromosome, ranges =
                                         IRanges(snvs$Start_Position,
                                                 snvs$End_Position), ref =
                                         snvs$Tumor_Seq_Allele1, alt =
                                         snvs$Tumor_Seq_Allele2)
  overlaps <- findOverlaps(vrange, genes)
  transcribed_variants <- rep("NA", nrow(dat))
  transcribed_variants[snv_index[queryHits(overlaps)]] <- as.character(decode(
    strand(genes[subjectHits(overlaps)])))

  #Match transcription and +, -, to account for reverse complement
  mut_type <- paste(dat$Tumor_Seq_Allele1, ">", dat$Tumor_Seq_Allele2,
                            sep = "")
  final_strand <- rep(NA, nrow(dat))
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- mut_type %in% forward_change
  final_strand[ind] <- ifelse(transcribed_variants[ind] == "+", "U", "T")
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- mut_type %in% rev_change
  final_strand[ind] <- ifelse(transcribed_variants[ind] == "-", "U", "T")

  dat[["Transcript_Strand"]] <- final_strand
  eval.parent(substitute(bay@variants <- dat))
  if (build_table) {
    dat_bagel <- methods::new("bagel", variants = drop_na_variants(
      dat, "Transcript_Strand"), count_tables = bay@count_tables,
      sample_annotations = bay@sample_annotations)
    tab <- build_custom_table(dat_bagel, variant_annotation =
                                  "Transcript_Strand", name =
                                  "Transcript_Strand", return_instead = FALSE)
    eval.parent(substitute(bay@count_tables <- tab))
  }
}

#' Add replication strand annotation to SNV variants based on bedgraph file
#'
#' @param bay Input bagel
#' @param rep_range A GRanges object with replication timing as metadata
#' @param build_table Automatically build a table from this annotation
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' annotate_replication_strand(bay, BAGEL::rep_range)
#' @export
annotate_replication_strand <- function(bay, rep_range, build_table = TRUE) {
  dat <- bay@variants
  snv_index <- which(dat$Variant_Type == "SNV")
  snvs <- subset_variant_by_type(dat, "SNV")

  #Create GRanges object to determine strand of variants within genes
  dat_range <- GenomicRanges::GRanges(seqnames = snvs$Chromosome, ranges =
                                         IRanges(snvs$Start_Position,
                                                 snvs$End_Position), ref =
                                         snvs$Tumor_Seq_Allele1, alt =
                                         snvs$Tumor_Seq_Allele2)
  overlaps <- GenomicRanges::findOverlaps(dat_range, rep_range)
  repl_variants <- rep("NA", length(dat_range))
  repl_variants[snv_index[S4Vectors::queryHits(overlaps)]] <-
    GenomicRanges::elementMetadata(BAGEL::rep_range)@listData[[1]][
      S4Vectors::subjectHits(overlaps)]
  dat[["Replication_Strand"]] <- repl_variants
  eval.parent(substitute(bay@variants <- dat))
  if (build_table) {
    dat_bagel <- methods::new("bagel", variants = drop_na_variants(
      dat, "Replication_Strand"), count_tables = bay@count_tables,
      sample_annotations = bay@sample_annotations)
    tab <- build_custom_table(dat_bagel, variant_annotation =
                                "Replication_Strand", name =
                                "Replication_Strand", data_factor =
                                factor(c("leading", "lagging")),
                              return_instead = FALSE)
    eval.parent(substitute(bay@count_tables <- tab))
  }
}
