#' @importFrom ggplot2 ggplot aes geom_bar theme theme_set geom_point aes_string
#' @importFrom ggplot2 facet_grid theme_bw xlab ylab element_blank element_text
#' @importFrom rlang .data
NULL

#' Return sample from bagel object
#'
#' @param bay Bagel object containing samples
#' @param sample_name Sample name to plot counts
#' @return Generates sample plot {no return}
#' @examples
#' bay <- readRDS(system.file("testdata", "bagel.rds", package = "BAGEL"))
#' plot_sample_counts(bay, get_sample_names(bay)[1])
#' @export
plot_sample_counts <- function(bay, sample_name) {
  if (length(sample_name) != 1) {
    stop("`please specify exactly one sample")
  }
  sample_number <- which(get_sample_names(bay = bay) == sample_name)
  sample <- bay@counts_table[, sample_number]
  plot_full(sample)
}

#' Complete Plotting of input data
#'
#' @param sample Single sample DataFrame
#' @return Generates sample plot {no return}
plot_full <- function(sample) {
  mut_summary <- table_96(sample)
  major <- table(mut_summary[, "Type"])
  df <- data.frame(major)
  mut <- data.frame(table(mut_summary$mutation))

  plot_major <- ggplot(df, aes_string(x = "Var1", y = "Freq",
                                               fill = "Var1")) + geom_bar(stat =
                                                               "identity") +
    theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) +
    xlab("Major Motif") + ylab("Counts") + theme(legend.position = "none")

  df2 <- data.frame(major = substr(mut$Var1, 1, 3), minor =
                      substr(mut$Var1, 5, 7), num = mut$Freq)
  df2$major <- as.character(df2$major)
  df2$minor <- as.character(df2$minor)
  plot_iris2 <- df2 %>% ggplot(aes_string(y = "num", x = "minor",
                                          fill = "major")) +
    ggplot2::facet_grid(~ major, scales = "free_x") +
    geom_bar(stat = "identity") + cowplot::background_grid(major = "y",
                                                           minor = "none") +
    cowplot::panel_border() + theme(axis.text.x = element_text(angle = 90,
                                                               vjust = 1,
                                                               hjust = 1)) +
    xlab("Minor Motif") + ylab("Counts") + ggplot2::labs(fill = "SNP")
  legend <- cowplot::get_legend(plot_iris2)
  plot_iris2 <- plot_iris2 + theme(legend.position = "none")

  theme_set(cowplot::theme_cowplot(font_size = 8)) # reduce default font size

  p <- cowplot::plot_grid(plot_iris2, legend, plot_major, labels = "AUTO",
                          ncol = 2, rel_widths = c(10, 1))
  return(p)
}

#' Plotting Signature Motif counts/spectra
#'
#' @param result S4 Result Object
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' plot_signatures(result)
#' @export
plot_signatures <- function(result) {
  signatures <- result@signatures
  groups <- reshape2::colsplit(rownames(signatures), "_", names = c("mutation",
                                                                    "context"))
  sig_names <- colnames(signatures)
  mutation <- groups[, "mutation"]

  signatures %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "Motif") %>%
    tidyr::gather("var", "val", -"Motif") %>%
    cbind(mutation) -> plot_dat

  plot_dat %>%
    ggplot(aes_string(y = "val", x = "Motif", fill = "mutation")) +
    geom_bar(stat = "identity") +
    facet_grid(factor(var, levels = unique(var), ordered = TRUE) ~ .,
               labeller = ggplot2::as_labeller(structure(sig_names, names =
                                                           unique(
                                                             plot_dat$var)))) +
    theme_bw() + xlab("Motifs") + ylab("Proportion") + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      strip.text.y = element_text(size = 7)) + ggplot2::scale_y_continuous(
        expand = c(0, 0)) -> p
  return(p)
}

#' Plotting Signature Motif counts/spectra
#'
#' @param result S4 Result Object
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' plot_signatures(result)
#' @export
plot_sample_reconstruction_error <- function(result, sample_number) {
  signatures <- result@bagel@counts_table[, sample_number, drop = FALSE]
  reconstructed <- reconstruct_sample(result, sample_number)
  error <- signatures - reconstructed
  colnames(error) <- "Error"

  groups <- reshape2::colsplit(rownames(signatures), "_", names = c("mutation",
                                                                    "context"))

  signatures %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "Motif") %>%
    tidyr::gather("var", "val", -"Motif") %>%
    cbind(groups[, "mutation"]) -> plot_dat

  reconstructed %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "Motif") %>%
    tidyr::gather("var", "val", -"Motif") %>%
    cbind(groups[, "mutation"]) -> reconstructed_dat

  error %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "Motif") %>%
    tidyr::gather("var", "val", -"Motif") %>%
    cbind(groups[, "mutation"]) -> error_dat

  plot_dat <- rbind(plot_dat, reconstructed_dat, error_dat)
  colnames(plot_dat) <- c("Motif", "var", "val", "mutation")

  plot_dat %>%
    ggplot(aes_string(y = "val", x = "Motif", fill = "mutation")) +
    geom_bar(stat = "identity") + theme_bw() + xlab("Motifs") +
    ylab("Proportion") + theme(axis.text.x = element_blank(),
                                                axis.ticks.x = element_blank(),
                                                strip.text.y =
                                                  element_text(size = 7)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::facet_wrap(~ var, drop = TRUE, scales = "fixed")
}

#' Plots signature weights for each sample
#'
#' @param result S4 Result Object
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @param label_samples Whether to display sample names or not
#' @param samples_plotted A subset of samples to plot
#' @param sort_samples Defaults to numerical, may be total 'counts', or a
#' specific signatures name (e.g. 'Signature1)
#' @param num_samples Number of sorted samples to plot
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' plot_exposures(result)
#' @export
plot_exposures <- function(result, proportional = TRUE, label_samples = TRUE,
                           samples_plotted = colnames(result@samples),
                           sort_samples = "numerical",
                           num_samples = length(colnames(result@samples))) {
  samples <- result@samples
  if (!all(samples_plotted %in% colnames(samples))) {
    stop(strwrap(prefix = " ", initial = "", "Some samples specified are
    not in this dataset, available samples are as follows: "), "\n",
         paste0(colnames(samples), collapse = "\n"))
  }
  samples <- samples[, samples_plotted]
  y_label <- "counts"
  if (proportional) {
    y_label <- "%"
    samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  }
  samples %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "make") %>%
    tidyr::gather("var", "val", -"make") -> plot_dat
  if (length(sort_samples == 1) && sort_samples == "numerical") {
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(gtools::mixedsort(plot_dat$var)))
  }else if (length(sort_samples == 1) && sort_samples == "counts") {
    counts <- colSums(samples)
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(plot_dat$var)
                           [order(counts, decreasing = TRUE)])
  }else {
    if (!all(sort_samples %in% rownames(samples))) {
      stop("Signature is not present in this result, please choose from: \n",
           paste0(rownames(samples), collapse = "\n"))
    }
    if (length(sort_samples) == 1) {
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[order(samples[
                                 sort_samples, ], decreasing = TRUE)])
    }
    else{
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[do.call(order,
                                                            c(decreasing = TRUE,
                                               as.data.frame(
                                                 t(samples[sort_samples, ]))))])
    }
  }
  samples_to_use <- levels(plot_dat$var)[seq_len(num_samples)]
  plot_dat <- plot_dat[which(plot_dat$var %in% samples_to_use), ]


  #Testing
  plot_dat_table <- plot_dat %>% dplyr::arrange(.data$make)
  if (all(!sort_samples %in% c("numerical", "counts"))) {
    plot_levels <- c(setdiff(unique(plot_dat_table$make), sort_samples),
                     rev(sort_samples))
    plot_dat_table$make <- factor(plot_dat_table$make, levels = plot_levels)
  }

  plot_dat_table %>%
    ggplot() + geom_bar(aes_string(y = "val", x = "var", fill = "make"),
                        stat = "identity") + theme_bw() + theme(legend.title =
                                                  element_blank(), axis.text.x =
                         element_text(angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.major.x = element_blank()) +
    xlab("Samples") + ylab(paste("Signatures (", y_label, ")", sep = "")) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) -> p
  if (!label_samples) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x =
                     element_blank())
  }
  return(p)
}

#' Plots signature weights for each sample gridded by single annotation
#'
#' @param result S4 Result Object
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @param label_samples Whether to display sample names or not
#' @param sort_samples Defaults to numerical, may be total 'counts', or a
#' specific signatures name (e.g. 'Signature1)
#' @param annotation Annotation to use for facet_wrap
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
#' plot_exposures(result)
#' @export
plot_exposures_by_annotation <- function(result, proportional = TRUE,
                                         label_samples = FALSE,
                                         sort_samples = "numerical",
                                         annotation, by_group = TRUE) {
  samples <- result@samples

  y_label <- "counts"
  if (proportional) {
    y_label <- "%"
    samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  }
  samples %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "make") %>%
    tidyr::gather("var", "val", -"make") -> plot_dat
  annot_dt <- result@bagel@sample_annotations
  plot_dat$annotation <- annot_dt[[annotation]][match(plot_dat$var,
                                                      annot_dt[["Samples"]])]
  #Normal alphabetical sorting of samples
  if (length(sort_samples == 1) && sort_samples == "numerical") {
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(gtools::mixedsort(plot_dat$var)))
    #Sorting of samples by counts
  }else if (length(sort_samples == 1) && sort_samples == "counts") {
    counts <- colSums(samples)
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(plot_dat$var)
                           [order(counts, decreasing = TRUE)])
    #Sorting of samples by counts of specific signature(s)
  }else {
    if (!all(sort_samples %in% rownames(samples))) {
      stop("Signature is not present in this result, please choose from: \n",
           paste0(rownames(samples), collapse = "\n"))
    }
    if (length(sort_samples) == 1) {
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[order(samples[
                                 sort_samples, ], decreasing = TRUE)])
    }
    else {
      plot_dat$var <- factor(plot_dat$var, levels = unique(plot_dat$var)[
        do.call(order, c(decreasing = TRUE, as.data.frame(t(samples[
          sort_samples, ]))))])
    }
  }
  plot_dat_table <- plot_dat %>% dplyr::arrange(.data$make)

  #If used, sorts indicated signature to the bottom so it's easier to see
  if (all(!sort_samples %in% c("numerical", "counts"))) {
    plot_levels <- c(setdiff(unique(plot_dat_table$make), sort_samples),
                     rev(sort_samples))
    plot_dat_table$make <- factor(plot_dat_table$make, levels = plot_levels)
  }
  if (by_group) {
    plot_dat_table %>%
      ggplot() + geom_bar(aes_string(y = "val", x = "var", fill = "make"),
                          stat = "identity") + theme_bw() + theme(
                            legend.title = element_blank(), axis.text.x =
              element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major.x = element_blank()) + xlab("Samples") +
      ylab(paste("Signatures (", y_label, ")", sep = "")) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) -> p
    p <- p + ggplot2::facet_wrap(~ annotation, drop = TRUE, scales = "free")
  } else {
    plot_dat_table %>%
      ggplot(aes_string(x = "annotation", y = "val", fill = "annotation")) +
      ggplot2::geom_boxplot() + theme_bw() + theme(legend.title =
                                                     element_blank(),
                                                   axis.text.x =
              element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major.x = element_blank()) + xlab("Annotation") +
      ylab(paste("Signatures (", y_label, ")", sep = "")) -> p
    p <- p + ggplot2::facet_wrap(~ make, drop = TRUE, scales = "free")
  }
  if (!label_samples) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x =
                     element_blank())
  }
  return(p)
}
