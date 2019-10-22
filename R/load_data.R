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
select_genome <- function(hg){
  #Choose genome build version
  #Keep for now as a helper function
  if (hg == 19){
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }else if (hg == 38){
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else{
    stop("That genome build is not currently supported")
    }
  return(g)
}

#' Loads a list of vcfs and converts to a combined data.table
#'
#' @param vcfs A list of vcf file locations
#' @return Returns a data.table of variants from vcfs
#' @examples
#' melanoma_vcf <- list.files(system.file("testdata", package = "BAGEL"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE)
#' melanoma <- BAGEL::vcfs_to_dt(vcfs = melanoma_vcf)
#' @export
vcfs_to_dt <- function(vcfs){
  vcf_list <- vector("list", length(vcfs))
  for (i in seq_len(length(vcfs))){
    vcf_list[[i]] <- BAGEL::vcf_to_dt(vcfs[i])
  }
  combined <- do.call("rbind", vcf_list)
  return(combined)
}

#' Loads a vcf and converts to data.table
#'
#' @param vcf_file Location of vcf file
#' @return Returns a data.table of variants from a vcf
#' @examples
#' luad1_vcf <- system.file("testdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "BAGEL")
#' luad1 <- BAGEL::vcf_to_dt(vcf = luad1_vcf)
#' @export
vcf_to_dt <- function(vcf_file){
  required_columns <- c("Chromosome", "Start_Position", "End_Position",
                       "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                       "Tumor_Sample_Barcode")

  vcf <- VariantAnnotation::readVcf(vcf_file)
  vcf_name <- basename(vcf_file)

  #Remove MultiAllelic Sites for now
  num_alleles <- lengths(VariantAnnotation::fixed(vcf)[, "ALT"])
  multi_allelic <- which(num_alleles != 1)
  if (length(multi_allelic > 0)){
    vcf <- vcf[-multi_allelic, ]
  }
  rows <- SummarizedExperiment::rowRanges(vcf)

  #Remove filtered rows
  pass <- rows$FILTER
  if (length(which(pass == "PASS")) < 1){
    warning(paste("No variants passed filtering, please review VCF file: ",
                  vcf_name, sep = ""))
  }
  rows <- rows[pass == "PASS", ]

  rows$REF <- vapply(as.character(rows$REF), substr, 1, 1, FUN.VALUE =
                       character(1))
  rows$ALT <- vapply(as.character(unlist(rows$ALT)), substr, 1, 1, FUN.VALUE =
                       character(1))

  dt <- cbind(data.table::as.data.table(rows),
              Tumor_Sample_Barcode = vcf_name) %>% dplyr::rename(
                "Tumor_Seq_Allele1" = "REF", "Tumor_Seq_Allele2" = "ALT",
                    "Chromosome" = "seqnames", "Start_Position" = "start",
                    "End_Position" = "end")
  variant_type <- rep(NA, nrow(dt))
  variant_type[which( (nchar(dt$Tumor_Seq_Allele1) == 1) &
                        (nchar(dt$Tumor_Seq_Allele2) == 1))] <- "SNP"
  dt <- dt[which(variant_type == "SNP"), ]

  if (!all(required_columns %in% names(dt))){
    warning("Some required columns missing")
  }
  dt <- dt[, required_columns]
  return(dt)
}

#' Loads a maf file and extracts the data.table of variants
#'
#' @param maf_file Location of maf file
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file=system.file("testdata", "TCGA_test.maf", package = "BAGEL")
#' maf = BAGEL::maf_to_dt(maf_file = maf_file)
#' @export
maf_to_dt <- function(maf_file){
  required_columns <- c("Chromosome", "Start_Position", "End_Position",
                       "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                       "Tumor_Sample_Barcode")
  maf <- maftools::read.maf(maf_file)
  dat <- maf@data
  if (!all(required_columns %in% names(dat))){
    warning("Some required columns missing")
  }
  dat <- dat[which(dat$Variant_Type == "SNP"), ]
  if (nrow(dat) == 0){
    warning(paste("No variants found in maf file:", basename(maf_file)))
  }
  dat <- dat[, required_columns]
  return(dat)
}
