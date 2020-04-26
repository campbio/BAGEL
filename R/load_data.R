# #' @importFrom SummarizedExperiment rowRanges
#' @importFrom methods is
NULL

#' Helper function to load HG19 or HG38
#'
#' @param hg Which human genome build number to return
#' @return Returns BSgenome of given version
#' @examples
#' g <- select_genome(hg = 38)
#'
#' @export
select_genome <- function(hg) {
  #Choose genome build version
  #Keep for now as a helper function
  if (hg == 19) {
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }else if (hg == 38) {
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else{
    stop("That genome build is not currently supported")
    }
  return(g)
}

#' Chooses the correct function to load in input based on class or file
#' extension
#'
#' @param input VCF, MAF, file.vcf, file.maf, or list of inputs
#' @param name Optional name for better error-reporting
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @param auto_fix_errors Attempt to automatically fix file formatting errors
#' @param verbose Show file list progress
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad_vcf_file <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad <- BAGEL::auto_to_bagel_dt(input = luad_vcf_file)
#'
#' luad_vcf_file <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad_vcf <- VariantAnnotation::readVcf(luad_vcf_file)
#' luad_vcf_name <- basename(luad_vcf_file)
#' luad <- BAGEL::auto_to_bagel_dt(input = luad_vcf, name = luad_vcf_name)
#'
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = maftools::read.maf(maf_file)
#' dt = BAGEL::maf_to_dt(maf)
#' maf_dt = BAGEL::auto_to_bagel_dt(input = dt, filter = FALSE)
#'
#' melanoma_vcfs <- list.files(system.file("testdata", package = "BAGEL"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE)
#' melanoma <- BAGEL::auto_to_bagel_dt(input = melanoma_vcfs)
#' @export
auto_to_bagel_dt <- function(input, name = NULL, filter = TRUE, only_snp = TRUE,
                             extra_fields = NULL, auto_fix_errors = TRUE,
                             verbose = TRUE) {
  if (length(input) > 1 && is(input, "vector")) {
    if(!is(input, "list")) {
      input <- as.list(input)
    }
    input_list <- vector("list", length(input))
    pb <- utils::txtProgressBar(min = 0, max = length(input_list), initial = 0,
                                style = 3)
    for (i in seq_len(length(input))) {
      utils::setTxtProgressBar(pb, i, )
      if (verbose) {
        print(paste("Sample number: ", i, "; Sample name: ", input[[i]],
                    sep = ""))
      }
      input_list[[i]] <- BAGEL::auto_to_bagel_dt(input = input[[i]],
                                                 name = name, filter = filter,
                                                 only_snp = only_snp,
                                                 extra_fields = extra_fields,
                                                 verbose = verbose)
    }
    dt <- do.call("rbind", input_list)
  } else if (is(input, "CollapsedVCF")) {
    dt <- vcf_to_dt(vcf = input, vcf_name = name, filter = filter,
                    only_snp = only_snp, extra_fields = extra_fields)
  } else if (is(input, "MAF")) {
    dt <- maf_to_dt(maf = input, maf_name = name, filter = filter,
                    only_snp = only_snp, extra_fields = extra_fields)
  } else if (is(input, "data.frame")) {
    dt <- dt_to_bagel_dt(dt = input, dt_name = name, filter = filter,
                         only_snp = only_snp, extra_fields = extra_fields)
  } else if (is(input, "character")) {
    if (tools::file_ext(input) == "vcf") {
      dt <- vcf_file_to_dt(vcf_file = input, filter = filter,
                           only_snp = only_snp, extra_fields = extra_fields,
                           auto_fix_errors = TRUE)
    } else if (tools::file_ext(input) == "maf") {
      dt <- maf_file_to_dt(maf_file = input, filter = filter,
                           only_snp = only_snp, extra_fields = extra_fields)
    } else if (!tools::file_ext(input) %in% c("vcf", "maf")) {
      stop(paste("Input: ", input, " could not be parsed.", sep = ""))
    }
  } else {
    stop(paste("There is no function to read this input. Input must be VCF, ",
    "MAF, file.vcf, file.maf, data.frame or a list of inputs. ",
    "Instead it is a ", class(input), sep = ""))
  }
  return(dt)
}

#' Loads a list of vcfs and converts to a combined data.table
#'
#' @param vcf_files A list of vcf file locations
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @param verbose Show file list progress
#' @return Returns a data.table of variants from vcfs
#' @examples
#' melanoma_vcfs <- list.files(system.file("testdata", package = "BAGEL"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE)
#' melanoma <- BAGEL::vcf_files_to_dt(vcf_files = melanoma_vcfs)
#' @export
vcf_files_to_dt <- function(vcf_files, filter = TRUE, only_snp = TRUE,
                       extra_fields = NULL, verbose = FALSE) {
  vcf_list <- vector("list", length(vcf_files))
  pb <- utils::txtProgressBar(min = 0, max = length(vcf_list), initial = 0,
                              style = 3)
  for (i in seq_len(length(vcf_files))) {
    utils::setTxtProgressBar(pb, i, )
    if (verbose) {
      print(paste("Sample number: ", i, "; Sample name: ", vcf_files[i],
                sep = ""))
    }
    vcf_list[[i]] <- BAGEL::vcf_file_to_dt(vcf_file = vcf_files[i],
                                           filter = filter, only_snp = only_snp,
                                           extra_fields = extra_fields)
  }
  combined <- do.call("rbind", vcf_list)
  return(combined)
}

#' Converts a loaded vcf object to data.table
#'
#' @param vcf Location of vcf file
#' @param vcf_name Name of the sample or vcf
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad_vcf_file <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad_vcf <- VariantAnnotation::readVcf(luad_vcf_file)
#' luad_vcf_name <- basename(luad_vcf_file)
#' luad <- BAGEL::vcf_to_dt(vcf = luad_vcf, vcf_name = luad_vcf_name)
#' @export
vcf_to_dt <- function(vcf, vcf_name = NULL, filter = TRUE, only_snp = TRUE,
                      extra_fields = NULL) {
  used_fields <- c(used_fields(), extra_fields)

  #Remove MultiAllelic Sites for now
  num_alleles <- lengths(VariantAnnotation::fixed(vcf)[, "ALT"])
  multi_allelic <- which(num_alleles != 1)
  if (length(multi_allelic) > 0) {
    vcf <- vcf[-multi_allelic, ]

  }

  rows <- SummarizedExperiment::rowRanges(vcf)

  if (filter) {
    #Remove filtered rows
    pass <- which(rows$FILTER == "PASS")
    if (length(pass) == 0) {
      warning(paste("No variants passed filtering, please review VCF
                           file: ", vcf_name, sep = ""))
      return(NULL)
    }
    rows <- rows[pass, ]
  }

  rows$REF <- vapply(as.character(rows$REF), substr, 1, 1,
                          FUN.VALUE = character(1))
  rows$ALT <- vapply(as.character(unlist(rows$ALT)), substr, 1, 1,
                          FUN.VALUE = character(1))

  if ("INFO" %in% used_fields) {
    info_field <- VariantAnnotation::info(vcf)
    info_field <- info_field[pass, , drop = FALSE]
    dt <- cbind(data.table::as.data.table(rows),
                Tumor_Sample_Barcode = vcf_name,
                data.table::as.data.table(info_field)) %>% dplyr::rename(
                  "Tumor_Seq_Allele1" = "REF", "Tumor_Seq_Allele2" = "ALT",
                      "Chromosome" = "seqnames", "Start_Position" = "start",
                      "End_Position" = "end")
  } else {
    dt <- cbind(data.table::as.data.table(rows),
                Tumor_Sample_Barcode = vcf_name) %>% dplyr::rename(
                  "Tumor_Seq_Allele1" = "REF", "Tumor_Seq_Allele2" = "ALT",
                  "Chromosome" = "seqnames", "Start_Position" = "start",
                  "End_Position" = "end")
  }

  dt <- add_variant_type(dt)

  if (only_snp) {
    variant_type <- rep(NA, nrow(dt))
    variant_type[which((nchar(dt$Tumor_Seq_Allele1) == 1) &
                          (nchar(dt$Tumor_Seq_Allele2) == 1))] <- "SNP"
    dt <- dt[which(variant_type == "SNP"), ]
  }

  if (length(used_fields) > 1 || used_fields[1] != FALSE) {
    if (!all(used_fields %in% names(dt))) {
      warning("Some required columns missing")
    }
    dt <- dt[, used_fields, with = FALSE]
  }

  #For some reason non-variants are included (e.g. T>T), remove them
  non_variant <- which(dt$Tumor_Seq_Allele1 == dt$Tumor_Seq_Allele2)
  if (length(non_variant) > 0) {
    dt <- dt[-non_variant, ]
  }

  #Drop factor levels which cause problems down the line
  dt[["Chromosome"]] <- as.character(dt[["Chromosome"]])
  GenomeInfoDb::seqlevelsStyle(dt$Chromosome) <- "UCSC"
  return(dt)
}

#' Loads a vcf and converts to data.table
#'
#' @param vcf_file Location of vcf file
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @param auto_fix_errors Attempt to automatically fix file formatting errors
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad_vcf <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad <- BAGEL::vcf_file_to_dt(vcf_file = luad_vcf)
#' @export
vcf_file_to_dt <- function(vcf_file, filter = TRUE, only_snp = TRUE,
                           extra_fields = NULL, auto_fix_errors = TRUE) {

  vcf <- try(VariantAnnotation::readVcf(vcf_file), silent = TRUE)
  vcf_name <- basename(vcf_file)

  if (class(vcf) == "try-error" && auto_fix_errors) {
    alt_input <- utils::read.table(vcf_file, stringsAsFactors = FALSE,
                            check.names = FALSE, comment.char = "", skip = 7,
                            header = TRUE)
    sample_header_name <- names(alt_input[10])
    vcf_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT")
    if (!all(vcf_columns %in% names(alt_input))) {
      stop(paste("VCF File: ", vcf_file,
                 " could not be recovered, please review.",
                 " \nAdditional information: \n", vcf[1], sep = ""))
    }
    non_letters <- rawToChar(as.raw(c(32:47, 58:64, 91, 93:96, 123:126)),
                             multiple = TRUE)
    malformed_rows <- c(which(alt_input$REF == 0), which(alt_input$ALT == 0),
                        grep("#", alt_input$`#CHROM`))
    if (length(malformed_rows) > 0) {
      alt_input <- alt_input[-malformed_rows, ]
    } else {
      stop(paste("VCF File: ", vcf_file,
                 " could not be recovered, please review.",
                 " \nAdditional information: \n", vcf[1], sep = ""))
    }
    alt_input[, "End_Position"] <- alt_input[, "POS"]
    alt_input <- alt_input[, colnames(alt_input)[c(1:2, ncol(alt_input),
                                                   4:(ncol(alt_input) - 1))]]
    dt <- cbind(data.table::as.data.table(alt_input[, c("#CHROM", "POS",
                                                        "End_Position", "REF",
                                                        "ALT", "QUAL",
                                                        "FILTER")]),
                Tumor_Sample_Barcode = vcf_name, data.table::as.data.table(
                  alt_input[, c("INFO", "FORMAT", sample_header_name), drop =
                              FALSE])) %>% dplyr::rename(
                                "Tumor_Seq_Allele1" = "REF",
                                "Tumor_Seq_Allele2" = "ALT", "Chromosome" =
                                  "#CHROM", "Start_Position" = "POS",
                                "Alleles" = sample_header_name)
    warning(paste("\nVCF File: ", vcf_file,
               " is malformed but could be recovered, review optional.",
               " \nAdditional information: \n", vcf[1],
               sep = ""))
    return(dt_to_bagel_dt(dt = dt, dt_name = vcf_name, filter = filter,
                          only_snp = only_snp, extra_fields = extra_fields))
  } else if (class(vcf) == "try-error") {
    stop(paste("VCF File: ", vcf_file,
                  " is malformed and auto-recovery is disabled, please review.",
                  " \nAdditional information: \n", vcf[1],
                  sep = ""))
  } else {
    return(vcf_to_dt(vcf = vcf, filter = filter, vcf_name = vcf_name,
                     only_snp = only_snp, extra_fields = extra_fields))
  }
}

#' Extracts the data.table of variants from a maftools maf object
#'
#' @param maf Maf object loaded by maftools::read.maf()
#' @param maf_name Name of the sample or maf
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = maftools::read.maf(maf_file)
#' maf_dt = BAGEL::maf_to_dt(maf = maf, maf_name = "test", filter = FALSE)
#' @export
maf_to_dt <- function(maf, maf_name = NULL, filter = TRUE, only_snp =
                              TRUE, extra_fields = NULL) {
  dt <- rbind(maf@data, maf@maf.silent)
  return(dt_to_bagel_dt(dt = dt, dt_name = maf_name, filter = filter,
                        only_snp = only_snp, extra_fields = extra_fields))
}

#' Converts a data.table to a filtered BAGEL data.table
#'
#' @param dt data.table input by user or vcf or maf functions to be filtered
#' @param dt_name Name of the vcf or maf
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = maftools::read.maf(maf_file)
#' dt = BAGEL::maf_to_dt(maf)
#' maf_dt = BAGEL::dt_to_bagel_dt(dt = dt, filter = FALSE)
#' @export
dt_to_bagel_dt <- function(dt, dt_name = NULL, filter = TRUE, only_snp = TRUE,
                           extra_fields = NULL) {
  used_fields <- c(used_fields(), extra_fields)
  if (is(dt, "data.frame") && !is(dt, "data.table")) {
    warning(paste(file_type, ": ", dt_name,
                  " is a data.frame but not a data.table, automatically fixing",
                  sep = ""))
    dt <- data.table::as.data.table(dt)
    file_type <- "data.table"
  } else if (is(dt, "data.table")) {
    file_type <- "data.table"
  } else {
    stop(paste(deparse(substitute(dt)), "/", dt_name,
               ": needs to be a data.frame or data.table it is a", class(dt),
               sep = ""))
  }
  if (!all(tolower(used_fields) %in% tolower(colnames(dt)))) {
    stop(paste("Required column(s) ",
                  used_fields[which(!tolower(used_fields) %in%
                                      tolower(colnames(dt)))], " missing in ",
                  file_type, ": ", dt_name, sep = ""))
  }
  for (i in seq_along(colnames(dt))) {
    dt_col <- colnames(dt)[i]
    if (any(tolower(dt_col) == tolower(used_fields)) && !dt_col %in%
        used_fields) {
      colnames(dt)[i] <- used_fields[which(grepl(dt_col, used_fields,
                                                 ignore.case = TRUE))]
      warning(paste("Column ",
                 dt_col, " had the wrong case and was automatically fixed ",
                 file_type, ": ", dt_name, sep = ""))
    }
  }


  if (filter) {
    pass <- which(dt$FILTER == "PASS")
    if (length(pass) == 0) {
      warning("No variants passed filtering")
    }
    dt <- dt[pass, ]
  }

  if (only_snp) {
    if (!is.null(dt$Variant_Type)) {
      snp_vars <- which(dt$Variant_Type == "SNP")
      if (length(snp_vars) > 0) {
        data.table::set(dt, snp_vars, "Variant_Type", "SNV")
      }
      dt <- dt[which(dt$Variant_Type == "SNV"), ]
    } else {
      dt <- add_variant_type(dt)
    }
    if (nrow(dt) == 0) {
      warning(paste("No variants found in ", file_type, ":", dt_name,
                    sep = ""))
    }
  } else {
    dt <- add_variant_type(dt)
  }
  dt <- dt[, used_fields, with = FALSE]

  #For some reason non-variants are included (e.g. T>T), remove them
  non_variant <- which(dt$Tumor_Seq_Allele1 == dt$Tumor_Seq_Allele2)
  if (length(non_variant) > 0) {
    dt <- dt[-non_variant, ]
  }
  #Drop factor levels which cause problems down the line
  dt[["Chromosome"]] <- as.character(dt[["Chromosome"]])
  GenomeInfoDb::seqlevelsStyle(dt$Chromosome) <- "UCSC"

  return(dt)
}

#' Loads a maf file and extracts the data.table of variants
#'
#' @param maf_file Location of maf file
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param extra_fields Which additional fields to extract
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = BAGEL::maf_file_to_dt(maf_file = maf_file)
#' @export
maf_file_to_dt <- function(maf_file, filter = TRUE, only_snp = TRUE,
                           extra_fields = NULL) {
  maf <- maftools::read.maf(maf_file)
  maf_name <- basename(maf_file)
  return(maf_to_dt(maf = maf, maf_name = maf_name, filter = filter,
                   only_snp = only_snp, extra_fields = extra_fields))
}

used_fields <- function() {
  return(c("Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Variant_Type"))
}
