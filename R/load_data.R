#' @importFrom SummarizedExperiment seqnames rowRanges
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

#' Loads a list of vcfs and converts to a combined data.table
#'
#' @param vcf_files A list of vcf file locations
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param used_fields Which fields to extract
#' @return Returns a data.table of variants from vcfs
#' @examples
#' melanoma_vcfs <- list.files(system.file("testdata", package = "BAGEL"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE)
#' melanoma <- BAGEL::vcf_files_to_dt(vcf_files = melanoma_vcfs)
#' @export
vcf_files_to_dt <- function(vcf_files, filter = TRUE, only_snp = TRUE,
                       used_fields = c("Chromosome", "Start_Position",
                                       "End_Position",
                                       "Tumor_Seq_Allele1",
                                       "Tumor_Seq_Allele2",
                                       "Tumor_Sample_Barcode")) {
  vcf_list <- vector("list", length(vcf_files))
  pb <- utils::txtProgressBar(min = 0, max = length(vcf_list), initial = 0,
                              style = 3)
  for (i in seq_len(length(vcf_files))) {
    utils::setTxtProgressBar(pb, i, )
    vcf_list[[i]] <- BAGEL::vcf_file_to_dt(vcf_file = vcf_files[i],
                                           filter = filter, only_snp = only_snp,
                                           used_fields = used_fields)
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
#' @param used_fields Which fields to extract
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad_vcf_file <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad_vcf <- VariantAnnotation::readVcf(luad_vcf_file)
#' luad_vcf_name <- basename(luad_vcf_file)
#' luad <- BAGEL::vcf_to_dt(vcf = luad_vcf, vcf_name = luad_vcf_name)
#' @export
vcf_to_dt <- function(vcf, vcf_name, filter = TRUE, only_snp = TRUE,
                      used_fields = c("Chromosome", "Start_Position",
                                            "End_Position", "Tumor_Seq_Allele1",
                                            "Tumor_Seq_Allele2",
                                            "Tumor_Sample_Barcode")) {

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

  rows$REF <- vapply(as.character(rows$REF), substr, 1, 1, FUN.VALUE =
                       character(1))
  rows$ALT <- vapply(as.character(unlist(rows$ALT)), substr, 1, 1, FUN.VALUE =
                       character(1))

  if ("INFO" %in% used_fields) {
    info_field <- VariantAnnotation::info(vcf)
    info_field <- info_field[-multi_allelic, , drop = FALSE]
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
  return(dt)
}

#' Loads a vcf and converts to data.table
#'
#' @param vcf_file Location of vcf file
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param used_fields Which fields to extract
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad_vcf <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad <- BAGEL::vcf_file_to_dt(vcf_file = luad_vcf)
#' @export
vcf_file_to_dt <- function(vcf_file, filter = TRUE, only_snp = TRUE,
                           used_fields = c("Chromosome", "Start_Position",
                                                "End_Position",
                                                "Tumor_Seq_Allele1",
                                                "Tumor_Seq_Allele2",
                                                "Tumor_Sample_Barcode")) {
  vcf <- VariantAnnotation::readVcf(vcf_file)
  vcf_name <- basename(vcf_file)
  return(vcf_to_dt(vcf = vcf, filter = filter, vcf_name = vcf_name,
            only_snp = only_snp, used_fields = used_fields))
}

#' Extracts the data.table of variants from a maftools maf object
#'
#' @param maf Maf object loaded by maftools::read.maf()
#' @param maf_name Name of the sample or maf
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param used_fields Which fields to extract
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = maftools::read.maf(maf_file)
#' maf_dt = BAGEL::maf_to_dt(maf = maf)
#' @export
maf_to_dt <- function(maf, maf_name = NULL, filter = TRUE, only_snp = TRUE,
                           used_fields = c("Chromosome", "Start_Position",
                                                "End_Position",
                                                "Tumor_Seq_Allele1",
                                                "Tumor_Seq_Allele2",
                                                "Tumor_Sample_Barcode")) {
  dt <- rbind(maf@data, maf@maf.silent)

  if (filter) {
    pass <- which(dt$FILTER == "PASS")
    if (length(pass) == 0) {
      warning("No variants passed filtering")
    }
    dt <- dt[pass, ]
  }

  if (only_snp) {
    dt <- dt[which(dt$Variant_Type == "SNP"), ]
    if (nrow(dt) == 0) {
      warning(paste("No variants found in maf file:", maf_name))
    }
  }

  if (length(used_fields) > 1 || used_fields[1] != FALSE) {
    if (!all(used_fields %in% names(dt))) {
      warning("Some required columns missing")
    }
    dt <- dt[, used_fields, with = FALSE]
  }
  return(dt)
}

#' Loads a maf file and extracts the data.table of variants
#'
#' @param maf_file Location of maf file
#' @param filter Filter to only passed variants
#' @param only_snp Filter only non-snp variants
#' @param used_fields Which fields to extract
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "public_TCGA.LUSC.maf", package = "BAGEL")
#' maf = BAGEL::maf_file_to_dt(maf_file = maf_file)
#' @export
maf_file_to_dt <- function(maf_file, filter = TRUE, only_snp = TRUE,
                           used_fields = c("Chromosome", "Start_Position",
                                                "End_Position",
                                                "Tumor_Seq_Allele1",
                                                "Tumor_Seq_Allele2",
                                                "Tumor_Sample_Barcode")) {
  maf <- maftools::read.maf(maf_file)
  maf_name <- basename(maf_file)
  return(maf_to_dt(maf = maf, maf_name = maf_name, filter = filter,
                   only_snp = only_snp, used_fields = used_fields))
}
