#' @importFrom NMF nmf
NULL

#' Discovers signatures and weights from a table of counts using NMF
#'
#' @param input A bagel object or counts table
#' @param table_name Name of table used for signature discovery
#' @param num_signatures Number of signatures to discover, k
#' @param method Discovery of new signatures using either LDA or NMF
#' @param seed Seed for reproducible signature discovery
#' @param nstart Number of independent runs with optimal chosen (lda only)
#' @param par_cores Number of parallel cores to use (NMF only)
#' @return Returns a result object with results and input object (if bagel)
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' build_standard_table(bay, "SNV96")
#' discover_signatures(input = bay, table_name = "SNV96",
#' num_signatures = 3, method = "nmf", seed = 12345, nstart = 1)
#' @export
discover_signatures <- function(input, table_name, num_signatures, method="lda",
                            seed = 1, nstart = 1, par_cores = FALSE) {
  if (!methods::is(input, "bagel")) {
    if (!methods::is(input, "matrix")) {
      stop("Input to discover_signatures must be a bagel object or a matrix")
    }
    bagel <- methods::new("bagel")
    bagel@count_tables@table_list[[1]] <- input
    bagel@count_tables@table_name[1] <- table_name
    input <- bagel
  }
  counts_table <- extract_count_table(input, table_name)

  #Determine if samples are present and can be used to scale weights
  used_samples <- which(input@variants$sample %in%
                          colnames(counts_table))
  if (length(used_samples) == 0) {
    warning(strwrap(prefix = " ", initial = "", "No samples overlap with
                      counts table, exposures will not be scaled by sample
                      counts."))
  } else {
    sample_counts <- table(input@variants$sample[used_samples])
    matched <- match(colnames(counts_table), names(sample_counts))
  }
  if (method == "lda") {
    counts_table <- t(counts_table)
    if (is.null(seed)) {
      control <- list(nstart = nstart)
    } else {
      control <- list(seed = (seq_len(nstart) - 1) + seed, nstart = nstart)
    }
    lda_out <- topicmodels::LDA(counts_table, num_signatures, control = control)
    lda_sigs <- exp(t(lda_out@beta))
    rownames(lda_sigs) <- colnames(counts_table)
    colnames(lda_sigs) <- paste("Signature", seq_len(num_signatures), sep = "")

    weights <- t(lda_out@gamma)
    rownames(weights) <- paste("Signature", seq_len(num_signatures), sep = "")
    colnames(weights) <- rownames(counts_table)

    # Multiply Weights by sample counts
    if (length(used_samples) != 0) {
      weights <- sweep(weights, 2, sample_counts[matched], FUN = "*")
    }
    lda_result <- methods::new("Result", signatures = lda_sigs,
                               exposures = weights, type = "LDA", bagel = input,
                               log_lik = stats::median(lda_out@loglikelihood),
                               perplexity = topicmodels::perplexity(lda_out))
    return(lda_result)
  } else if (method == "nmf") {
    #Needed to prevent error with entirely zero rows
    epsilon <- 0.00000001
    if (par_cores) {
      decomp <- NMF::nmf(counts_table + epsilon, num_signatures, seed = seed,
                         nrun = nstart, .options = paste("p", par_cores,
                                                         sep = ""))
    } else {
      decomp <- NMF::nmf(counts_table + epsilon, num_signatures, seed = seed,
                         nrun = nstart)
    }
    rownames(decomp@fit@H) <- paste("Signature", seq_len(num_signatures),
                                 sep = "")
    colnames(decomp@fit@W) <- paste("Signature", seq_len(num_signatures),
                                    sep = "")
    nmf_result <- methods::new("Result", signatures = decomp@fit@W,
                               exposures = decomp@fit@H, type = "NMF",
                               bagel = input, log_lik = decomp@residuals)
    nmf_result@signatures <- sweep(nmf_result@signatures, 2,
                                   colSums(nmf_result@signatures), FUN = "/")
    nmf_result@exposures <- sweep(nmf_result@exposures, 2,
                                   colSums(nmf_result@exposures), FUN = "/")

    # Multiply Weights by sample counts
    if (length(used_samples) != 0) {
      nmf_result@exposures <- sweep(nmf_result@exposures, 2,
                                  sample_counts[matched], FUN = "*")
    }
    return(nmf_result)
  } else{
    stop("That method is not supported. Use lda or nmf to generate signatures.")
  }
}

kld <- function(a, b) {
  return(sum(a * log2(a / b)))
}

#' Compare two vectors similarity based on Jensen-Shannon Divergence
#'
#' @param p First vector
#' @param q Second vector
#' @return Returns Jensen-Shannon Divergence of the two vectors
#' @examples
#' p <- c(0.2, 0.3, 0.4, 0.5, 0.6)
#' q <- c(0.6, 0.5, 0.4, 0.3, 0.2)
#' jsd(p, q)
#' @export
jsd <- function(p, q) {
  epsilon <- 0.0000001
  p <- p + epsilon
  q <- q + epsilon
  m <- (p + q) / 2
  jsd <- 0.5 * kld(p, m) + 0.5 * kld(q, m)
  return(jsd)
}

sig_compare <- function(sig1, sig2, threshold=0.9) {
  sig1_names <- colnames(sig1)
  sig2_names <- colnames(sig2)
  if (nrow(sig1) != nrow(sig2)) {
    stop("Signatures must have the same motifs")
  }
  matches <- matrix(nrow = ncol(sig1), ncol = ncol(sig2))
  for (i in seq_len(ncol(sig1))) {
    for (j in seq_len(ncol(sig2))) {
      matches[i, j] <- 1 - jsd(sig1[, i], sig2[, j])
    }
  }
  comparison <- NULL
  for (row in seq_len(nrow(matches))) {
    line <- which(matches[row, ] > threshold)
    if (length(line) > 0) {
      for (match in line) {
        comparison <- rbind(comparison, c(matches[row, match], row, match,
                                          sig1_names[row], sig2_names[match]))
      }
    }
  }
  if (is.null(comparison)) {
    stop("No matches found, try lowering threshold.")
  }
  comparison <- data.frame(comparison, stringsAsFactors = FALSE)
  colnames(comparison) <- c("cor", "xindex", "yindex", "xcol", "ycol")
  comparison$cor <- as.numeric(comparison$cor)
  comparison$xindex <- as.numeric(comparison$xindex)
  comparison$yindex <- as.numeric(comparison$yindex)
  return(comparison)
}

#' Compare two result files to find similar signatures
#'
#' @param result Result to compare
#' @param other_result Second result
#' @param threshold threshold for similarity
#' @param result_name title for plot of first result signatures
#' @param other_result_name title for plot of second result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' compare_results(res, res, threshold = 0.8)
#' @export
compare_results <- function(result, other_result,
                            threshold = 0.9, result_name =
                              deparse(substitute(result)), other_result_name =
                              deparse(substitute(other_result))) {
  signatures <- result@signatures
  comparison <- sig_compare(signatures, other_result@signatures, threshold)
  result_subset <- methods::new("Result",
                      signatures = result@signatures[, comparison$xindex,
                                                     drop = FALSE], exposures =
                        matrix(), type = "NMF")
  other_subset <- methods::new("Result",
                      signatures = other_result@signatures[, comparison$yindex,
                                                            drop = FALSE],
                      exposures = matrix(), type = "NMF")
  result_plot <- BAGEL::plot_signatures(result_subset)
  result_plot <- result_plot + ggplot2::ggtitle(result_name)
  cosmic_plot <- BAGEL::plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle(other_result_name)
  gridExtra::grid.arrange(result_plot, cosmic_plot, ncol = 2)
  return(comparison)
}

#' Compare a result object to COSMIC V3 Signatures; Select exome or genome for
#' SNV and only genome for DBS or Indel classes
#'
#' @param result Result to compare
#' @param variant_class Compare to SNV, DBS, or Indel
#' @param sample_type exome (SNV only) or genome
#' @param threshold threshold for similarity
#' @param result_name title for plot user result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' compare_cosmic_v3(res, "SNV", "genome", threshold = 0.8)
#' @export
compare_cosmic_v3 <- function(result, variant_class, sample_type,
                              threshold = 0.9, result_name =
                              deparse(substitute(result))) {
  if (sample_type == "exome") {
    if (variant_class %in% c("snv", "SNV", "SNV96", "SBS")) {
      cosmic_res <- cosmic_v3_snv_sigs_exome
    } else {
      stop(paste("Only SNV class is available for whole-exome, please choose",
                 " `genome` for DBS or Indel", sep = ""))
    }
  } else if (sample_type == "genome") {
    if (variant_class %in% c("snv", "SNV", "SNV96", "SBS")) {
      cosmic_res <- cosmic_v3_snv_sigs
    } else if (variant_class %in% c("DBS", "dbs", "doublet")) {
      cosmic_res <- cosmic_v3_dbs_sigs
    } else if (variant_class %in% c("INDEL", "Indel", "indel", "ind", "IND",
                                  "ID")) {
      cosmic_res <- cosmic_v3_indel_sigs
    } else {
      stop("Only SNV, DBS, and Indel classes are supported")
    }
  } else {
    stop("Sample type must be exome or genome")
  }
  signatures <- result@signatures
  comparison <- sig_compare(signatures, cosmic_res@signatures, threshold)
  result_subset <- methods::new(
    "Result", signatures = result@signatures[,
                                             comparison$xindex, drop = FALSE],
    exposures = matrix(), type = "NMF")
  other_subset <- methods::new(
    "Result", signatures = cosmic_res@signatures[, comparison$yindex,
                                                 drop = FALSE],
    exposures = matrix(), type = "NMF")
  result_plot <- BAGEL::plot_signatures(result_subset)
  result_plot <- result_plot + ggplot2::ggtitle(result_name)
  cosmic_plot <- BAGEL::plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle(paste("COSMIC Signatures v3",
                                                      variant_class, " ",
                                                      sample_type, sep = ""))
  gridExtra::grid.arrange(result_plot, cosmic_plot, ncol = 2)
  return(comparison)
}

#' Compare a result object to COSMIC V2 SNV Signatures (combination whole-exome
#' and whole-genome)
#'
#' @param result Result to compare
#' @param threshold threshold for similarity
#' @param result_name title for plot user result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' compare_cosmic_v2(res, threshold = 0.8)
#' @export
compare_cosmic_v2 <- function(result, threshold = 0.9, result_name =
                                deparse(substitute(result))) {
  signatures <- result@signatures
  comparison <- sig_compare(signatures, cosmic_v2_sigs@signatures, threshold)
  result_subset <- methods::new("Result",
                                signatures =
                                  result@signatures[, comparison$xindex, drop =
                                                      FALSE], exposures =
                                  matrix(), type = "NMF")
  other_subset <- methods::new("Result",
                               signatures =
                                 cosmic_v2_sigs@signatures[, comparison$yindex,
                                                           drop = FALSE],
                               exposures = matrix(), type = "NMF")
  result_plot <- BAGEL::plot_signatures(result_subset)
  legend <- cowplot::get_legend(result_plot)
  result_plot <- result_plot + ggplot2::ggtitle(result_name) +
    theme(legend.position = "none")
  cosmic_plot <- BAGEL::plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle("COSMIC Signatures v2") +
    theme(legend.position = "none")
  gridExtra::grid.arrange(result_plot, cosmic_plot, legend, ncol = 3,
                          widths = c(0.4, 0.4, 0.2))
  return(comparison)
}

#' Input a cancer subtype to return a list of related COSMIC signatures
#'
#' @param tumor_type Cancer subtype to view related signatures
#' @return Returns signatures related to a partial string match
#' @examples cosmic_v2_subtype_map ("lung")
#' @export
cosmic_v2_subtype_map <- function(tumor_type) {
  subtypes <- c("adrenocortical carcinoma", "all", "aml", "bladder", "breast",
               "cervix", "chondrosarcoma", "cll", "colorectum", "glioblastoma",
               "glioma low grade", "head and neck", "kidney chromophobe",
               "kidney clear cell", "kidney papillary", "liver", "lung adeno",
               "lung  small cell", "lung  squamous", "lymphoma b-cell",
               "lymphoma hodgkin", "medulloblastoma", "melanoma", "myeloma",
               "nasopharyngeal carcinoma", "neuroblastoma", "oesophagus",
               "oral gingivo-buccal squamous", "osteosarcoma", "ovary",
               "pancreas", "paraganglioma", "pilocytic astrocytoma", "prostate",
               "stomach", "thyroid", "urothelial carcinoma", "uterine carcinoma"
               , "uterine carcinosarcoma", "uveal melanoma")
  present_sig <- list(
    c(1, 2, 4, 5, 6, 13, 18), c(1, 2, 5, 13), c(1, 5), c(1, 2, 5, 10, 13),
    c(1, 2, 3, 5, 6, 8, 10, 13, 17, 18, 20, 26, 30), c(1, 2, 5, 6, 10, 13, 26),
    c(1, 5), c(1, 2, 5, 9, 13), c(1, 5, 6, 10), c(1, 5, 11), c(1, 5, 6, 14),
    c(1, 2, 4, 5, 7, 13), c(1, 5, 6), c(1, 5, 6, 27), c(1, 2, 5, 13),
    c(1, 4, 5, 6, 12, 16, 17, 22, 23, 24), c(1, 2, 4, 5, 6, 13, 17),
    c(1, 4, 5, 15), c(1, 2, 4, 5, 13), c(1, 2, 5, 9, 13, 17), c(1, 5, 25),
    c(1, 5, 8), c(1, 5, 7, 11, 17), c(1, 2, 5, 13), c(1, 2, 5, 6, 13),
    c(1, 5, 18), c(1, 2, 4, 5, 6, 13, 17), c(1, 2, 5, 7, 13, 29),
    c(1, 2, 5, 6, 13, 30), c(1, 5), c(1, 2, 3, 5, 6, 13), c(1, 5), c(1, 5, 19),
    c(1, 5, 6), c(1, 2, 5, 13, 15, 17, 18, 20, 21, 26, 28), c(1, 2, 5, 13),
    c(1, 2, 5, 13, 22), c(1, 2, 5, 6, 10, 13, 14, 26), c(1, 2, 5, 6, 10, 13),
    c(1, 5, 6)
  )
  partial <- grep(tumor_type, subtypes)
  for (i in seq_len(length(partial))) {
    print(subtypes[partial[i]])
    print(present_sig[[partial[i]]])
  }
}

#' LDA prediction of samples based on existing signatures
#'
#' @param bagel Input samples to predit signature weights
#' @param table_name Name of table used for posterior prediction.
#' Must match the table type used to generate the prediction signatures
#' @param signature_res Signatures to use for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", decompTumor2Sig, and deconstructSigs
#' @param signatures_to_use Which signatures in set to use (default all)
#' @param seed Seed to use for reproducible results, set to null to disable
#' @param verbose Whether to show intermediate results
#' @return Results a result object containing signatures and sample weights
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' build_standard_table(bay, "SNV96")
#' predict_exposure(bay, "SNV96", BAGEL::cosmic_v2_sigs, algorithm = "lda")
#' @export
predict_exposure <- function(bagel, table_name, signature_res, algorithm,
                             signatures_to_use = seq_len(ncol(
                               signature_res@signatures)), seed = 1,
                             verbose = FALSE) {
  signature <- signature_res@signatures[, signatures_to_use]
  counts_table <- extract_count_table(bagel, table_name)
  present_samples <- which(colSums(counts_table) > 0)
  counts_table <- counts_table[, present_samples]

  if (algorithm %in% c("lda_posterior", "lda", "lda_post")) {
    ifelse(is.null(seed),
           lda_res <- lda_posterior(counts_table = counts_table,
                                    signature = signature, max.iter = 100,
                                    verbose = verbose),
           lda_res <- withr::with_seed(1, lda_posterior(
                        counts_table = counts_table, signature = signature,
                        max.iter = 100, verbose = verbose)))
    exposures <- t(lda_res$samp_sig_prob_mat)
    type_name <- "posterior_LDA"
  }else if (algorithm %in% c("decomp", "decompTumor2Sig")) {
    ifelse(is.null(seed),
           decomp_res <- predict_decompTumor2Sig(counts_table, signature),
           decomp_res <- withr::with_seed(seed, predict_decompTumor2Sig(
                           counts_table, signature)))
    exposures <- t(do.call(rbind, decomp_res))
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
    type_name <- "decompTumor2Sig"
  }else if (algorithm %in% c("ds", "deconstruct", "deconstructSigs")) {
    ifelse(is.null(seed),
           sigs.input <- deconstructSigs:: mut.to.sigs.input(mut.ref =
                                                               bagel@variants,
                                      sample.id = "Tumor_Sample_Barcode",
                                      chr = "Chromosome",
                                      pos = "Start_Position",
                                      ref = "Tumor_Seq_Allele1",
                                      alt = "Tumor_Seq_Allele2",
                                      bsg = bagel@genome),
           sigs.input <- withr::with_seed(
             deconstructSigs::mut.to.sigs.input(mut.ref = bagel@variants,
                                             sample.id = "Tumor_Sample_Barcode",
                                             chr = "Chromosome",
                                             pos = "Start_Position",
                                             ref = "Tumor_Seq_Allele1",
                                             alt = "Tumor_Seq_Allele2",
                                             bsg = bagel@genome)))
    sig_all <- t(signature)
    middle <- unlist(lapply(strsplit(colnames(sig_all), "_"), "[", 1))
    context <- lapply(strsplit(colnames(sig_all), "_"), "[", 2)
    first <- unlist(lapply(context, substr, 1, 1))
    last <- unlist(lapply(context, substr, 3, 3))
    new_cols <- paste(first, "[", middle, "]", last, sep="")
    colnames(sig_all) <- new_cols

    ds_res <- sapply(rownames(sigs.input), function(x) {
      ds_result <- whichSignatures(tumor.ref = sigs.input,
                                   contexts.needed = TRUE,
                                   signatures.limit = ncol(signature),
                                   tri.counts.method = "default",
                                   sample.id = x, signatures.ref = sig_all)
      return(as.matrix(ds_result$weights))
    exposures <- ds_res
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
    type <- "deconstructSigs"
    })
  } else {
    stop("Type must be lda or decomp")
  }
  result <- methods::new("Result", signatures = signature,
                                       exposures = exposures,
                                       type = type_name, bagel = bagel)

  # Multiply Weights by sample counts
  used_samples <- which(bagel@variants$sample %in%
                          colnames(counts_table))
  if (length(used_samples) == 0) {
    warning(strwrap(prefix = " ", initial = "", "No samples overlap with
                      counts table, exposures will not be scaled by sample
                      counts."))
  } else {
    sample_counts <- table(bagel@variants$sample[used_samples])
    matched <- match(colnames(counts_table), names(sample_counts))
    result@exposures <- sweep(result@exposures, 2, sample_counts[matched],
                              FUN = "*")
  }
  return(result)
}

lda_posterior <- function(counts_table, signature, max.iter = 100,
                          theta = 0.1, verbose) {
  k <- ncol(signature) # number of signatures/topics
  num_samples <- ncol(counts_table) # number of samples

  if (length(theta) == 1) {
    theta <- rep(theta, k) # symmetric singular value converted to vector
  }
  sample_count_sums <- colSums(counts_table)

  # Initialize signature proportion matrix
  samp_sig_prob_mat <- matrix(NA, nrow = num_samples, ncol = k)
  sig_mut_counts <- matrix(NA, nrow = num_samples, ncol = k)
  rownames(samp_sig_prob_mat) <-
    rownames(sig_mut_counts) <- colnames(counts_table)
  colnames(samp_sig_prob_mat) <-
    colnames(sig_mut_counts) <- colnames(signature)

  for (s in seq_len(num_samples)) {
    sig_mut_counts[s, ] <- base::tabulate(sample(x = seq_len(k), size =
                                                   sample_count_sums[s],
                                                 replace = TRUE), k)
  }

  # Update signature proportion matrix
  if (verbose) {
    print("Calculating Signature Proportions")
  }
  for (i in seq_len(max.iter)) {
    for (s in seq_len(num_samples)) {
      #updating each mutation probability to reassign to a signature
      log_prob_mut_reassignment <-
        digamma(sig_mut_counts[s, ] + theta) -
        digamma(sample_count_sums[s] + sum(theta))
      #updating present sample topic probability
      sig_sample_weights <- t(signature + 1e-20) *
        exp(log_prob_mut_reassignment) # avoid 0 in norm
      sig_sample_weights <- sweep(sig_sample_weights, MARGIN = 2, STATS =
                                    colSums(sig_sample_weights), FUN = "/")
      #assigned counts for a topic for a sample
      updated_topic_motifs <- counts_table[, s] * t(sig_sample_weights)

      # Update nN.SbyT[s, ] sample counts assigned to signature
      sig_mut_counts[s, ] <- colSums(updated_topic_motifs)

      # Update p.SbyT[s, ]
      samp_sig_prob_mat[s, ] <- (sig_mut_counts[s, ]) / (sample_count_sums[s])
    }
    # Update theta
    theta <- MCMCprecision::fit_dirichlet(x = samp_sig_prob_mat)$alpha
    if (verbose) {
      print(theta)
    }
  }
  return(list(samp_sig_prob_mat = samp_sig_prob_mat, theta.poster = theta))
}

predict_decompTumor2Sig <- function(sample_mat, signature_mat){
  #Alexandrov-type prediction
  input_signatures_normalized <- apply(signature_mat, 2,
                                       function(x){x / sum(x)})
  signatures <- split(input_signatures_normalized,
                      col(input_signatures_normalized))
  signatures_ref <- decompTumor2Sig::readAlexandrovSignatures()
  ns <- as.matrix(row.names(signature_mat))
  ns <- apply(ns, 1, function(x){stringr::str_c(substr(x,5,5), '[',
                                                substr(x, 1, 3), ']',
                                                substr(x, 7, 7))})
  signatures <- lapply(signatures, setNames, ns)

  input_samples_normalized <- apply(sample_mat, 2, function(x){x / sum(x)})
  input_samples1 <- split(input_samples_normalized,
                          col(input_samples_normalized))
  genomes <- lapply(input_samples1, setNames, ns)

  sample_weight_mat <- decompTumor2Sig::decomposeTumorGenomes(genomes, signatures, verbose=FALSE)
  return(sample_weight_mat)
}

#placeholder
.multi_modal_discovery <- function(bay, num_signatures, motif96_name,
                                  rflank_name, lflank_name, max.iter=125,
                                  seed=123) {
  motif96 <- extract_count_table(bay, motif96_name)
  rflank <- extract_count_table(bay, rflank_name)
  lflank <- extract_count_table(bay, lflank_name)
  print(dim(motif96))
  print(dim(rflank))
  print(dim(lflank))
}

whichSignatures = function(tumor.ref = NA,
                           sample.id,
                           signatures.ref,
                           associated = c(),
                           signatures.limit = NA,
                           signature.cutoff = 0.06,
                           contexts.needed = FALSE,
                           tri.counts.method = "default") {
  if(class(tumor.ref) == 'matrix'){
    stop(paste('Input tumor.ref needs to be a data frame or location of input text file', sep = ''))
  }

  if(exists("tumor.ref", mode = "list")){
    tumor     <- tumor.ref
    if(contexts.needed == TRUE){
      tumor   <- deconstructSigs::getTriContextFraction(mut.counts.ref = tumor,
                                                        trimer.counts.method =
                                                          tri.counts.method)
    }
  } else {
    if(file.exists(tumor.ref)){
      tumor   <- utils::read.table(tumor.ref, sep = "\t", header = TRUE,
                                   as.is = TRUE, check.names = FALSE)
      if(contexts.needed == TRUE){
        tumor <- deconstructSigs::getTriContextFraction(tumor,
                                                        trimer.counts.method =
                                                          tri.counts.method)
      }
    } else {
      print("tumor.ref is neither a file nor a loaded data frame")
    }
  }

  if (missing(sample.id) && nrow(tumor) == 1) {
    sample.id = rownames(tumor)[1]
  }
  # Take patient id given
  tumor <- as.matrix(tumor)
  if(!sample.id %in% rownames(tumor)){
    stop(paste(sample.id, " not found in rownames of tumor.ref", sep = ''))
  }
  tumor <- subset(tumor, rownames(tumor) == sample.id)
  if(round(rowSums(tumor), digits = 1) != 1){
    stop(paste('Sample: ', sample.id, ' is not normalized\n', 'Consider using "contexts.needed = TRUE"', sep = ' '))
  }
  signatures <- signatures.ref

  signatures    <- as.matrix(signatures)
  original.sigs <- signatures

  # Check column names are formatted correctly
  if(length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))){
    colnames(tumor) <- deconstructSigs::changeColumnNames(colnames(tumor))
    if(length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))){
      stop("Check column names on input file")
    }
  }

  # Ensure that columns in tumor match the order of those in signatures
  tumor <- tumor[,colnames(signatures), drop = FALSE]

  #Take a subset of the signatures
  if(!is.null(associated)){
    signatures <- signatures[rownames(signatures) %in% associated, ]
  }

  if(is.na(signatures.limit)){
    signatures.limit <- nrow(signatures)
  }

  #Set the weights matrix to 0
  weights         <- matrix(0, nrow = nrow(tumor), ncol = nrow(signatures), dimnames = list(rownames(tumor), rownames(signatures)))

  seed            <- deconstructSigs::findSeed(tumor, signatures)
  weights[seed]   <- 1
  w               <- weights*10

  error_diff      <- Inf
  error_threshold <- 1e-3

  num <- 0
  while(error_diff > error_threshold){
    num        <- num + 1
    #print(num)
    error_pre  <- deconstructSigs::getError(tumor, signatures, w)
    if(error_pre == 0){
      break
    }
    w          <- deconstructSigs::updateW_GR(tumor, signatures, w,
                                              signatures.limit =
                                                signatures.limit)
    error_post <- deconstructSigs::getError(tumor, signatures, w)
    error_diff <- (error_pre-error_post)/error_pre
  }

  weights <- w/sum(w)
  unknown <- 0

  ## filtering on a given threshold value (0.06 default)
  weights[weights < signature.cutoff ] <- 0
  unknown <- 1 - sum(weights)

  product <- weights %*% signatures
  diff    <- tumor - product

  x       <- matrix(data = 0, nrow = 1, ncol = nrow(original.sigs),
                    dimnames = list(rownames(weights), rownames(original.sigs)))
  x       <- data.frame(x)
  x[colnames(weights)] <- weights
  weights <- x

  out        <- list(weights, tumor, product, diff, unknown)
  names(out) <- c("weights", "tumor", "product", "diff", "unknown")
  return(out)
}

#' Generate result_grid from bagel based on annotation and range of k
#'
#' @param bagel Input bagel to generate grid from
#' @param table_name Name of table used for signature discovery
#' @param discovery_type Algorithm for signature discovery
#' @param annotation Sample annotation to split results into
#' @param k_start Lower range of number of signatures for discovery
#' @param k_end Upper range of number of signatures for discovery
#' @param n_start Number of times to discover signatures and compare based on
#' posterior loglikihood
#' @param seed Seed to use for reproducible results, set to null to disable
#' @param par_cores Number of parallel cores to use (NMF only)
#' @param verbose Whether to output loop iterations
#' @return Results a result object containing signatures and sample weights
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel_snv96.rds", package = "BAGEL"))
#' grid <- generate_result_grid(bay, "SNV96", "nmf", k_start = 2, k_end = 5)
#' @export
generate_result_grid <- function(bagel, table_name, discovery_type = "lda",
                                 annotation = NA, k_start, k_end, n_start = 1,
                                 seed = NULL, par_cores = FALSE,
                                 verbose = FALSE) {
  result_grid <- methods::new("Result_Grid")

  #Set Parameters
  params <- data.table::data.table("discovery_type" = discovery_type,
                                   "annotation_used" = annotation, "k_start" =
                                     k_start, "k_end" = k_end,
                                   "total_num_samples" =
                                     nrow(bagel@sample_annotations),
                                   "nstart" = n_start, seed = seed)
  result_grid@grid_params <- params

  #Initialize grid_table and result_list
  grid_table <- data.table::data.table(annotation = character(), k =
                                         numeric(), num_samples = numeric(),
                                       reconstruction_error = numeric())
  result_list <- list()
  list_elem <- 1

  #Generate and set result_list
  if (!is.na(annotation)) {
    annot_samples <- bagel@sample_annotations$Samples
    annot <- bagel@sample_annotations$Tumor_Type
    annot_names <- unique(annot)
    num_annotation <- length(annot_names)
  } else {
    annot_names <- NA
    num_annotation <- 1
  }

  #Define new bagels
  for (i in 1:num_annotation) {
    if (!is.na(annotation)) {
      if (verbose) {
        cat(paste("Current Annotation: ", annot_names[i], "\n", sep = ""))
      }
      cur_ind <- which(annot == annot_names[i])
      cur_annot_samples <- annot_samples[cur_ind]
      cur_annot_variants <- bagel@variants[which(
        bagel@variants$sample %in% cur_annot_samples), ]

      cur_bagel <- methods::new("bagel", variants = cur_annot_variants,
                       sample_annotations =
                         bagel@sample_annotations[cur_ind, ],
                       count_tables = subset_count_tables(bagel,
                                                          cur_annot_samples))
    } else {
      cur_bagel <- bagel
      cur_annot_samples <- unique(bagel@variants$sample)
    }
    #Used for reconstruction error
    cur_counts <- extract_count_table(cur_bagel, table_name)

    #Define new results
    for (cur_k in k_start:k_end) {
      cur_result <- discover_signatures(input = cur_bagel, table_name =
                                          table_name, num_signatures = cur_k,
                                        method = discovery_type, nstart =
                                          n_start, seed = seed, par_cores =
                                          par_cores)
      result_list[[list_elem]] <- cur_result
      list_elem <- list_elem + 1

      recon_error <- mean(sapply(seq_len(ncol(cur_counts)), function(x)
        mean((cur_counts[, x, drop = FALSE] -
                reconstruct_sample(cur_result, x))^2))^2)

      grid_table <- rbind(grid_table, data.table::data.table(
        annotation = annot_names[i], k = cur_k, num_samples =
          length(cur_annot_samples), reconstruction_error = recon_error))
    }
  }
  result_grid@result_list <- result_list
  result_grid@grid_table <- grid_table
  return(result_grid)
}

reconstruct_sample <- function(result, sample_number) {
  reconstruction <- matrix(apply(sweep(result@signatures, 2,
                                       result@exposures[, sample_number,
                                                      drop = FALSE], FUN = "*"),
                                 1, sum), dimnames =
                             list(rownames(result@signatures), "Reconstructed"))
  return(reconstruction)
}

#' Automatic filtering of signatures for exposure prediction gridded across
#' specific annotation
#'
#' @param bagel Input samples to predit signature weights
#' @param table_name Name of table used for posterior prediction (e.g. SNV96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", decompTumor2Sig, and deconstructSigs
#' @param sample_annotation Annotation to grid across, if none given,
#' prediction subsetting on all samples together
#' @param min_exists Threshold to consider a signature active in a sample
#' @param proportion_samples Threshold of samples to consider a signature
#' active in the cohort
#' @param rare_exposure A sample will be considered active in the cohort if at
#' least one sample has more than this threshold proportion
#' @param verbose Print current annotation value being predicted on
#' @param combine_res Automatically combines a list of annotation results
#' into a single result object with zero exposure values for signatures not
#' found in a given annotation's set of samples
#' @param seed Seed to use for reproducible results, set to null to disable
#' @return Results a list of results, one per unique annotation value, if no
#' annotation value is given, returns a single result for all samples, or
#' combines into a single result if combines_res = TRUE
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel_annot.rds", package = "BAGEL"))
#' auto_predict_grid(bay, "SNV96", BAGEL::cosmic_v2_sigs, "lda",
#' "Tumor_Subtypes")
#'
#' auto_predict_grid(bay, "SNV96", BAGEL::cosmic_v2_sigs, "lda")
#' @export
auto_predict_grid <- function(bagel, table_name, signature_res, algorithm,
                              sample_annotation = NULL, min_exists = 0.05,
                              proportion_samples = 0.25, rare_exposure = 0.4,
                              verbose = TRUE, combine_res = TRUE, seed = 1) {
  if (is.null(sample_annotation)) {
    combine_res = FALSE
    result = auto_subset_sigs(bagel = bagel, table_name =
                       table_name, signature_res =
                       signature_res, algorithm = algorithm,
                       min_exists = min_exists, proportion_samples =
                       proportion_samples, rare_exposure =
                       rare_exposure, seed = seed)
  } else {
    available_annotations <- setdiff(colnames(bagel@sample_annotations),
                                     "Samples")
    if (!sample_annotation %in% available_annotations) {
      stop(paste0("Sample annotation ", sample_annotation, " not found, ",
                 "available annotations: ", available_annotations))
    }
    annot <- unique(bagel@sample_annotations[[sample_annotation]])
    result <- list()
    for (i in seq_along(annot)) {
      if (verbose) {
        print(as.character(annot[i]))
      }
      current_bagel <- BAGEL::subset_bagel_by_annotation(bagel, annot_col =
                                                           sample_annotation,
                                                         annot_name = annot[i])
      current_predicted <- auto_subset_sigs(bagel = current_bagel, table_name =
                                              table_name, signature_res =
                                              signature_res, min_exists =
                                              min_exists, proportion_samples =
                                              proportion_samples,
                                            rare_exposure = rare_exposure,
                                            algorithm = algorithm, seed = seed)
      result[[as.character(annot[i])]] <- current_predicted
    }
  }
  if (combine_res) {
    result = combine_predict_grid(result, bagel, signature_res)
  }
  return(result)
}

#' Automatic filtering of inactive signatures
#'
#' @param bagel Input samples to predit signature weights
#' @param table_name Name of table used for posterior prediction (e.g. SNV96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", decompTumor2Sig, and deconstructSigs
#' @param min_exists Threshold to consider a signature active in a sample
#' @param proportion_samples Threshold of samples to consider a signature
#' active in the cohort
#' @param rare_exposure A sample will be considered active in the cohort if at
#' least one sample has more than this threshold proportion
#' @param seed Seed to use for reproducible results, set to null to disable
#' @return Results a result object containing automatically subset signatures
#' and corresponding sample weights
auto_subset_sigs <- function(bagel, table_name, signature_res, algorithm,
                             min_exists = 0.05, proportion_samples = 0.25,
                             rare_exposure = 0.4, seed = 1) {
  test_predicted <- predict_exposure(bagel = bagel, table_name = table_name,
                                    signature_res = signature_res,
                                    algorithm = algorithm)
  exposures <- test_predicted@exposures
  num_samples <- ncol(exposures)
  exposures <- sweep(exposures, 2, colSums(exposures), "/")
  to_use <- as.numeric(which(apply(exposures, 1, function(x)
    sum(x > min_exists) / num_samples) > proportion_samples |
      apply(exposures, 1, max) > rare_exposure))
  final_inferred <- predict_exposure(bagel = bagel, table_name = table_name,
                                     signature_res = signature_res,
                                     signatures_to_use = to_use,
                                     algorithm = algorithm, seed = seed)
  return(final_inferred)
}

#' Combine prediction grid list into a result object. Exposure values are zero
#' for samples in an annotation where that signature was not predicted
#'
#' @param grid_list A list of result objects from the prediction grid to
#' combine into a single result
#' @param bagel Input samples to predit signature weights
#' @param signature_res Signatures to automatically subset from for prediction
#' @return A result object combining all samples and signatures from a
#' prediction grid. Samples have zero exposure value for signatures not found
#' in that annotation type.
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel_annot.rds", package = "BAGEL"))
#' grid <- auto_predict_grid(bay, "SNV96", BAGEL::cosmic_v2_sigs, "lda",
#' "Tumor_Subtypes", combine_res = FALSE)
#' combined <- combine_predict_grid(grid, bay, BAGEL::cosmic_v2_sigs)
#' plot_exposures_by_annotation(combined, "Tumor_Subtypes")
#' @export
combine_predict_grid <- function(grid_list, bagel, signature_res) {
  sig_names <- NULL
  for(i in 1:length(grid_list)) {
    sig_names <- c(sig_names, rownames(grid_list[[i]]@exposures))
  }
  sig_names <- unique(sig_names)
  sig_names <- sig_names[order(sig_names)]

  comb <- NULL
  for(i in 1:length(grid_list)) {
    samp <- grid_list[[i]]@exposures
    missing <- sig_names[!sig_names %in% rownames(samp)]
    missing_mat <- matrix(0, length(missing), ncol(samp))
    rownames(missing_mat) <- missing
    samp <- rbind(samp, missing_mat)
    samp <- samp[order(rownames(samp)), ]
    comb <- cbind(comb, samp)
  }
  grid_res <- new("Result", bagel = bagel, exposures = comb,
                  signatures = signature_res@signatures[, sig_names])
}
