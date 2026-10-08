# ---------------------------------------------------------------------------
# lib/plotting.R - one figure-saving path and one palette for every stage.
#
# Exports:
#   save_figure()             write a paired PDF and PNG, cairo where available
#   suz_theme()               shared ggplot theme
#   suz_direction_colours     UP / DOWN / NS
#   suz_rna_class_colours     mRNA / lncRNA / miRNA / circRNA
#   suz_phase_colours         Early / Sustained / Late / Reversal
#   suz_direction_levels      DOWN at the base, UP on top
#
# Figures are written at --dpi, default 600.
# ---------------------------------------------------------------------------

SUZ_DPI_DEFAULT <- 600

suz_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1,
                                               colour = "grey30"),
      legend.key       = ggplot2::element_blank()
    )
}

#' Save a ggplot as a PDF/PNG pair.
#'
#' @param plot    A ggplot (or anything ggsave accepts).
#' @param out_dir Output directory.
#' @param stem    Filename without extension; "_<dpi>dpi" is appended to the PNG.
save_figure <- function(plot, out_dir, stem, width = 8, height = 6,
                        dpi = SUZ_DPI_DEFAULT, png = TRUE, pdf = TRUE) {
  pdf_dev <- if (requireNamespace("grDevices", quietly = TRUE) &&
                 capabilities("cairo")) grDevices::cairo_pdf else NULL

  if (pdf) {
    p <- file.path(out_dir, paste0(stem, ".pdf"))
    ok <- tryCatch({
      if (!is.null(pdf_dev)) {
        ggplot2::ggsave(p, plot, width = width, height = height,
                        device = pdf_dev, limitsize = FALSE)
      } else {
        ggplot2::ggsave(p, plot, width = width, height = height,
                        device = "pdf", useDingbats = FALSE, limitsize = FALSE)
      }
      TRUE
    }, error = function(e) {
      warning("cairo_pdf failed for ", stem, " (", conditionMessage(e),
              "); falling back to pdf()", call. = FALSE)
      ggplot2::ggsave(p, plot, width = width, height = height,
                      device = "pdf", limitsize = FALSE)
      TRUE
    })
    if (ok) message("[fig] ", basename(p))
  }

  if (png) {
    p <- file.path(out_dir, paste0(stem, "_", dpi, "dpi.png"))
    ggplot2::ggsave(p, plot, width = width, height = height, dpi = dpi,
                    limitsize = FALSE)
    message("[fig] ", basename(p))
  }
  invisible(out_dir)
}

#' Shared direction palette. UP/DOWN colours are consistent across every panel
#' so Affinity assembly does not have to reconcile them.
suz_direction_colours <- c(
  UP        = "#B2182B",
  DOWN      = "#2166AC",
  NS        = "grey70",
  Unchanged = "grey70"
)

#' Shared RNA-class palette (Figure 2B, Figure 3A/B).
suz_rna_class_colours <- c(
  mRNA             = "#1B7837",
  lncRNA           = "#762A83",
  `annotated lncRNA` = "#9970AB",
  `novel lncRNA`   = "#C2A5CF",
  miRNA            = "#E08214",
  circRNA          = "#35978F",
  other            = "grey60"
)

#' Shared temporal-phase palette (Early / Sustained / Late).
suz_phase_colours <- c(
  Early     = "#FDB863",
  Sustained = "#5AAE61",
  Late      = "#8073AC",
  Reversal  = "#D6604D"
)

#' Stacked bar ordering used by Figure 3B and S3C: DOWN at the bottom,
#' UP at the top, white in-bar counts (36b_ and 39e_ both insist on this).
suz_direction_levels <- c("DOWN", "UP")
