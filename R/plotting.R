#' @importFrom ggplot2 ggplot aes geom_bar theme theme_set geom_point aes_string
#' @importFrom ggplot2 facet_grid theme_bw xlab ylab element_blank element_text
#' @importFrom rlang .data
NULL

#' Return sample from bagel object
#'
#' @param bay Bagel object containing samples
#' @param sample_name Sample name to plot counts
#' @return Generates sample plot {no return}
plot_sample <- function(bay, sample_name){
  if (length(sample_name) != 1){
    stop("`please specify exactly one sample")
  }
  sample_number <- which(sample_names(bay = bay) == sample_name)
  sample <- bay@counts_table[, sample_number]
  plot_full(sample)
}

#' Complete Plotting of input data
#'
#' @param sample Single sample DataFrame
#' @return Generates sample plot {no return}
plot_full <- function(sample){
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
#' result <- readRDS(system.file("testdata", "res.rds", package = "MotifSig"))
#' plot_signatures(result)
#' @export
plot_signatures <- function(result){
  signatures <- result@signatures
  groups <- reshape2::colsplit(rownames(signatures), "_", names = c("mutation",
                                                                    "context"))
  sig_names <- colnames(signatures)
  mutation <- groups[, "mutation"]

  signatures %>% as.data.frame %>% tibble::rownames_to_column(var = "Motif") %>%
    tidyr::gather("var", "val", -"Motif") %>% cbind(mutation) -> plot_dat

  plot_dat %>% ggplot(aes_string(y = "val", x = "Motif", fill = "mutation")) +
    geom_bar(stat = "identity") + facet_grid(factor(var, levels = unique(var),
                                                    ordered = TRUE) ~ .,
                                             labeller = ggplot2::as_labeller(
                                               structure(sig_names, names =
                                                           unique(plot_dat$var
                                                                  )))) +
    theme_bw() +
    xlab("Motifs") + ylab("Counts") + theme(axis.text.x = element_blank(),
                                            axis.ticks.x = element_blank(),
                                            strip.text.y =
                                              element_text(size = 7)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}

#' Plots signature weights for each sample
#'
#' @param result S4 Result Object
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "MotifSig"))
#' plot_samples(result)
#' @export
plot_samples <- function(result, proportional = TRUE){
  samples <- result@samples
  y_label <- "absolute"
  if (proportional){
    y_label <- "%"
    samples <- apply(samples, 2, function(x){
      x * 100 / sum(x, na.rm = TRUE)
    })
  }
  a <- samples %>% as.data.frame %>%
    tibble::rownames_to_column(var = "make") %>% tidyr::gather("var", "val",
                                                               -"make")
  a$var <- factor(a$var, levels = unique(gtools::mixedsort(a$var)))
  a %>% dplyr::arrange(.data$make) %>% ggplot() +
    geom_bar(aes_string(y = "val", x = "var", fill = "make"), stat =
               "identity") + theme_bw() + theme(legend.title =
                                                  element_blank(), axis.text.x =
                         element_text(angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.major.x = element_blank()) +
    xlab("Samples") + ylab(paste("Signatures (", y_label, ")", sep = "")) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}
