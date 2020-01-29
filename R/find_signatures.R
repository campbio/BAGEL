#' Discovers signatures and weights from a table of counts using NMF
#'
#' @param input A bagel object or counts table
#' @param num_signatures Number of signatures to discover, k
#' @param method Discovery of new signatures using either LDA or NMF
#' @return Returns a result object with results and input object (if bagel)
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' find_signatures(bay, num_signatures = 4)
#' @export
find_signatures <- function(input, num_signatures, method="lda", seed = NA,
                            nstart = 1) {
  if (methods::is(input, "bagel")) {
    has_tables(input)
  }else{
    if (!methods::is(input, "matrix")) {
      stop("Input to find_signatures must be a bagel object or a matrix")
    }
    bagel <- methods::new("bagel")
    bagel@counts_table <- input
    input <- bagel
  }

  #Determine if samples are present and can be used to scale weights
  used_samples <- which(input@variants$Tumor_Sample_Barcode %in%
                          colnames(input@counts_table))
  if (length(used_samples) == 0) {
    warning(strwrap(prefix = " ", initial = "", "No samples overlap with
                      counts table, exposures will not be scaled by sample
                      counts."))
  } else {
    sample_counts <- table(input@variants$Tumor_Sample_Barcode[used_samples])
    matched <- match(colnames(input@counts_table), names(sample_counts))
  }
  if (method == "lda") {
    counts_table <- t(input@counts_table)
    if (is.na(seed)) {
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
                               samples = weights, type = "LDA", bagel = input,
                               log_lik = median(lda_out@loglikelihood))
    return(lda_result)
  } else if (method == "nmf") {
    decomp <- NNLM::nnmf(input@counts_table, num_signatures, rel.tol = 1e-5,
                         max.iter = 10000L)
    rownames(decomp$H) <- paste("Signature", seq_len(num_signatures),
                                 sep = "")
    colnames(decomp$W) <- paste("Signature", seq_len(num_signatures), sep = "")
    nmf_result <- methods::new("Result", signatures = decomp$W,
                               samples = decomp$H, type = "NMF", bagel = input)
    nmf_result@signatures <- sweep(nmf_result@signatures, 2,
                                   colSums(nmf_result@signatures), FUN = "/")
    nmf_result@samples <- sweep(nmf_result@samples, 2,
                                   colSums(nmf_result@samples), FUN = "/")

    # Multiply Weights by sample counts
    if (length(used_samples) != 0) {
      nmf_result@samples <- sweep(nmf_result@samples, 2,
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

#' Compare two result files or input one to compare to COSMIC
#'
#' @param result Result to compare
#' @param other_result Second result, leave blank to use cosmic results
#' @param threshold threshold for similarity
#' @param result_name title for plot of first result signatures
#' @param other_result_name title for plot of second result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' compare_results(res, threshold = 0.8)
#' @export
compare_results <- function(result, other_result = cosmic_result,
                            threshold = 0.9, result_name = "User Signatures 1",
                            other_result_name = "User Signatures 2") {
  signatures <- result@signatures
  comparison <- sig_compare(signatures, other_result@signatures, threshold)
  result_subset <- methods::new("Result",
                      signatures = result@signatures[, comparison$xindex,
                                                     drop = FALSE], samples =
                        matrix(), type = "NMF")
  other_subset <- methods::new("Result",
                      signatures = other_result@signatures[, comparison$yindex,
                                                            drop = FALSE],
                      samples = matrix(), type = "NMF")
  if (identical(other_result, cosmic_result)) {
    result_name <- "User Signatures"
    other_result_name <- "COSMIC Signatures"
  }
  result_plot <- BAGEL::plot_signatures(result_subset)
  result_plot <- result_plot + ggplot2::ggtitle(result_name)
  cosmic_plot <- BAGEL::plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle(other_result_name)
  gridExtra::grid.arrange(result_plot, cosmic_plot, ncol = 2)
  return(comparison)
}

#' Input a cancer subtype to return a list of related COSMIC signatures
#'
#' @param tumor_type Cancer subtype to view related signatures
#' @return Returns signatures related to a partial string match
#' @examples what_cosmic30_sigs("lung")
#' @export
what_cosmic30_sigs <- function(tumor_type) {
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

#' LDA prediction of samples based on existing signatures (default COSMIC)
#'
#' @param bagel Input samples to predit signature weights
#' @param signatures Signatures to use for prediction (default COSMIC)
#' @param signatures_to_use Which signatures in set to use (default all)
#' @param verbose Whether to show intermediate results
#' @return Results a result object containing signatures and sample weights
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' infer_signatures(bay)
#' @export
infer_signatures <- function(bagel, signatures=cosmic_result@signatures,
                          signatures_to_use = seq_len(ncol(signatures)),
                          verbose = FALSE) {
  has_tables(bagel)

  # Load sample counts matrix
  counts_matrix <- bagel@counts_table

  # convert data structures
  sig_name <- colnames(signatures[, signatures_to_use])
  sig_props <- as.matrix(signatures[, colnames(signatures) %in% sig_name])
  samples_counts <- as.matrix(counts_matrix)

  est_sig_prop <- function(samples_counts, sig_props, max.iter = 100,
                           theta = 0.1) {
    k <- ncol(sig_props) # number of signatures/topics
    num_samples <- ncol(samples_counts) # number of samples

    if (length(theta) == 1) {
      theta <- rep(theta, k) # symmetric singular value converted to vector
    }
    sample_count_sums <- colSums(samples_counts)

    # Initialize signature proportion matrix
    samp_sig_prob_mat <- matrix(NA, nrow = num_samples, ncol = k)
    sig_mut_counts <- matrix(NA, nrow = num_samples, ncol = k)
    rownames(samp_sig_prob_mat) <-
      rownames(sig_mut_counts) <- colnames(samples_counts)
    colnames(samp_sig_prob_mat) <-
      colnames(sig_mut_counts) <- colnames(sig_props)

    for (s in seq_len(num_samples)) {
      sig_mut_counts[s, ] <- base::tabulate(sample(x = seq_len(k), size =
                                              sample_count_sums[s], replace =
                                              TRUE), k)
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
        sig_sample_weights <- t(sig_props + 1e-20) *
          exp(log_prob_mut_reassignment) # avoid 0 in norm
        sig_sample_weights <- sweep(sig_sample_weights, MARGIN = 2, STATS =
                                       colSums(sig_sample_weights), FUN = "/")
        #assigned counts for a topic for a sample
        updated_topic_motifs <- samples_counts[, s] * t(sig_sample_weights)

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
  res2 <- est_sig_prop(samples_counts = samples_counts, sig_props = sig_props,
                       max.iter = 100)
  lda_posterior_result <- methods::new("Result", signatures =
                               signatures[, signatures_to_use], samples =
                               t(res2$samp_sig_prob_mat), type =
                               "posterior_LDA", bagel = bagel)

  # Multiply Weights by sample counts
  used_samples <- which(bagel@variants$Tumor_Sample_Barcode %in%
                          colnames(bagel@counts_table))
  if (length(used_samples) == 0) {
    warning(strwrap(prefix = " ", initial = "", "No samples overlap with
                      counts table, exposures will not be scaled by sample
                      counts."))
  } else {
    sample_counts <- table(bagel@variants$Tumor_Sample_Barcode[used_samples])
    matched <- match(colnames(bagel@counts_table), names(sample_counts))
    lda_posterior_result@samples <- sweep(lda_posterior_result@samples, 2,
                                          sample_counts[matched], FUN = "*")
  }
  return(lda_posterior_result)
}

has_tables <- function(bagel) {
  if (!methods::is(bagel, "bagel")) {
    stop(strwrap(prefix = " ", initial = "", "The input object is not a
    bagel object, please use new('bagel') to create one."))
  }
  if (length(bagel@counts_table) == 0) {
    stop(strwrap(prefix = " ", initial = "", "The counts table is
    either missing or malformed, please run create_tables prior to this
                         function."))
  }
}

#' Generate result_grid from bagel based on annotation and range of k
#'
#' @param bagel Input bagel to generate grid from
#' @param discovery_type Algorithm for signature discovery
#' @param annotation Sample annotation to split results into
#' @param k_start Lower range of number of signatures for discovery
#' @param k_end Upper range of number of signatures for discovery
#' @param num_iter Number of times to discover signatures and compare based on
#' loglikihood
#' @param seed Give a seed to generate reproducible signatures
#' @param verbose Whether to output loop iterations
#' @return Results a result object containing signatures and sample weights
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' generate_result_grid(bay, k_start = 2, k_end = 5)
#' @export
generate_result_grid <- function(bagel, discovery_type = "lda", annotation,
                                 k_start, k_end, n_start = 1, seed = NA,
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
  annot_samples <- bagel@sample_annotations$Samples
  annot <- bagel@sample_annotations$Tumor_Type
  annot_names <- unique(annot)
  num_annotation <- length(annot_names)

  #Define new bagels
  for (i in 1:num_annotation) {
    if (verbose) {
      cat(paste("Current Annotation: ", annot_names[i], "\n", sep = ""))
    }
    cur_ind <- which(annot == annot_names[i])
    cur_annot_samples <- annot_samples[cur_ind]
    cur_annot_variants <- bagel@variants[which(
      bagel@variants$Tumor_Sample_Barcode %in% cur_annot_samples), ]

    cur_bagel <- new("bagel", variants = cur_annot_variants,
                     sample_annotations =
                       bagel@sample_annotations[cur_ind, ],
                     counts_table = bagel@counts_table[, cur_ind])

    #Used for reconstruction error
    cur_counts <- cur_bagel@counts_table

    #Define new results
    for (cur_k in k_start:k_end) {
      cur_result <- find_signatures(input = cur_bagel, num_signatures = cur_k,
                                    method = discovery_type, nstart = n_start,
                                    seed = seed)
      result_list[[list_elem]] <- cur_result
      list_elem <- list_elem + 1

      recon_error <- mean(sapply(seq_len(cur_counts), function(x)
        mean((cur_counts[, x, drop = FALSE] - reconstruct_sample(
          cur_result, x))^2))^2)
      #recon_error <- mean(sapply(1:ncol(cur_counts), function(x)
      #  mean((cur_counts[, x, drop = FALSE] - reconstruct_sample(
      #    cur_result, x))^2)/sum(cur_counts[, x, drop = FALSE]))^2)

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
                                       result@samples[, sample_number,
                                                      drop = FALSE], FUN = "*"),
                                 1, sum), dimnames =
                             list(rownames(result@signatures), "Reconstructed"))
  return(reconstruction)
}
