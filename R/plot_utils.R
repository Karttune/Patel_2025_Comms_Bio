########################################################################################################################################################################################################
############################################################# Function for plotting a browser snapshot from a region ###################################################################################
########################################################################################################################################################################################################

plot_single_region <- function(bigwigs, locus, gene_annotation = NULL, exons = F, highlighted_regs = NULL, tcga = F, title = NULL, ncols = 2, cores = 16) {
  require(tidyCoverage)
  require(rtracklayer)
  require(BiocParallel)

  cores <- MulticoreParam(cores)

  tracks <- BigWigFileList(bigwigs)
  names(tracks) <- names(bigwigs)

  ce <- CoverageExperiment(
    tracks,
    locus,
    width = round(width(locus) / 10) * 10, # Rounding to the nearest multiple of ten:
    # scale = TRUE,
    # center = TRUE,
    BPPARAM = cores
  ) |>
    coarsen(window = 50) |>
    expand()

  if (tcga == T) { # If input data is TCGA ATAC, replicates are combined
    ce <- ce |> mutate(
      rep = str_extract(track, "T1|T2"),
      track = str_extract(track, "TCGA-[:alnum:]{2}-[:alnum:]{4}-[:alnum:]{3}_[:alnum:]{4}")
    )
  }

  p <- ce |>
    ggplot(aes(x = coord, y = coverage)) +
    geom_col(alpha = .75) +
    facet_wrap(~track, ncol = ncols) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("#f05d5e", "#0f7173")) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 8)) +
    coord_cartesian(xlim = c(start(locus), end(locus))) +
    labs(title = title) +
    xlab("Coordinate (bp)") +
    ylab("Coverage")

  if (!is.null(gene_annotation)) {
    # Adding gene annotation track if the file is given:
    tabix_list <- Rsamtools::scanTabix(gene_annotation, param = locus)
    gtf_cnames <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    
    df_list <- lapply(tabix_list, function(x) {
      if (length(x) > 0) {
        if (length(x) == 1) {
          # Hack to make sure that it also works for data frames with only one row
          # Adds an empty row and then removes it
          result <- paste(paste(x, collapse = "\n"), "\n", sep = "")
          result <- readr::read_delim(result, delim = "\t", col_names = gtf_cnames)[1, ]
        } else {
          result <- paste(x, collapse = "\n")
          result <- readr::read_delim(result, delim = "\t", col_names = gtf_cnames)
        }
      } else {
        # Return NULL if the nothing is returned from tabix file
        result <- NULL
      }
      return(result)
    })
    
    gene_reg <- df_list[[1]] |>
      filter(feature == "gene") |>
      mutate(
        row_n = row_number(),
        gene_name = str_extract(attribute, '(?<=gene_id ")[^"]+'),
        hgnc_symbol = str_extract(attribute, '(?<=gene_name ")[^"]+')
      ) |>
      rowwise() |>
      mutate(
        gene_dir = if_else(strand == "+", "first", "last"),
        tss = ifelse(strand == "+", start, end),
        labelx = case_when(
          start < start(locus) & end < end(locus) ~ end + (width(locus) * 0.01), # If start coordinate of gene is SMALLER than the region, labeling the END coordinate if it is within the region, also offsetting the label by 1% of the total plotting region width
          start > start(locus) & end > end(locus) ~ start - (width(locus) * 0.01), # If start coordinate of gene is LARGER than the region, labeling the START coordinate if it is within the region
          start > start(locus) & end < end(locus) & strand == "+" ~ start - (width(locus) * 0.01), # If the whole gene is WITHIN the region, labeling the START, + strand
          start > start(locus) & end < end(locus) & strand == "-" ~ end + (width(locus) * 0.01), # If the whole gene is WITHIN the region, labeling the START
          start < start(locus) & end > end(locus) ~ start(locus) + ((end(locus) - start(locus)) / 2)
        ), # If the gene is larger than the whole region, labeling the CENTER
        just = case_when(
          start < start(locus) & end < end(locus) ~ "left", # If start coordinate of gene is SMALLER than the region, labeling the END coordinate if it is within the region
          start > start(locus) & end > end(locus) ~ "right", # If start coordinate of gene is LARGER than the region, labeling the START coordinate if it is within the region
          start > start(locus) & end < end(locus) & strand == "+" ~ "right", # If the whole gene is WITHIN the region, labeling the START
          start > start(locus) & end < end(locus) & strand == "-" ~ "left", # If the whole gene is WITHIN the region, labeling the START
          start < start(locus) & end > end(locus) ~ "center"
        ),
        lab = if_else(just != "center", hgnc_symbol, NA),
        # aoff = 0 - (max(ce$coverage) / 10) - row_n
        aoff = 0 - row_n * (max(ce$coverage) * 0.1) # To position the gene annotation correctly, taking a fraction of the maximum y axis coverage value and adding it as a negative multiplier
      )

    p <- p + geom_segment(inherit.aes = F, data = gene_reg, mapping = aes(x = if_else(strand == "+", start, end), xend = if_else(strand == "+", end, start), y = aoff, yend = aoff), arrow = arrow(length = unit(0.1, "cm")), size = .1) +
      geom_text(inherit.aes = F, data = gene_reg, aes(x = labelx, y = aoff, label = lab, hjust = just), size = 2, fontface = "italic")
  }

  if (exons == T) {
    exons <- df_list[[1]] |> 
      filter(feature == "exon") |> 
      mutate(
        gene_name = str_extract(attribute, '(?<=gene_id ")[^"]+'),
        hgnc_symbol = str_extract(attribute, '(?<=gene_name ")[^"]+')
      ) |> 
      rename(exon_start = start, exon_end = end) |> 
      dplyr::select(exon_start, exon_end, gene_name) |> 
      left_join(gene_reg, by = "gene_name")
    
    p <- p + geom_rect(inherit.aes = F, data = exons, aes(xmin = exon_start, xmax = exon_end, ymin = aoff - aoff * .1, ymax = aoff + aoff * .1))
  }

  if (!is.null(highlighted_regs)) {
    p <- p + geom_rect(inherit.aes = F, data = highlighted_regs |> filter_by_overlaps(locusus) |> as_tibble(), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), alpha = .2, fill = "#efcf82", color = "gray40", linewidth = .1, linetype = "dashed")
  }

  return(p)
}

########################################################################################################################################################################################################
######################################################## Function for plotting a metaplot of bigwig signals using tidyCoverage #########################################################################
########################################################################################################################################################################################################

plot_metaplot <- function(bigwigs, loci, bin_size = 25, reg_width = 1000, title = NULL, ncols = 2, nrows = 2) {
  # Requires ggsci
  require(ggsci)
  require(tidyCoverage)
  require(rtracklayer)
  require(parallel)
  require(BiocParallel)

  cores <- MulticoreParam(16)

  tracks <- BigWigFileList(bigwigs)

  if (!is.null(names(tracks))) { # If given bigwigs are names, naming the tracks too:
    names(tracks) <- names(bigwigs)
  } else {
    names(tracks) <- 1:length(bigwigs)
  }

  ce <- CoverageExperiment(tracks,
    loci,
    width = reg_width,
    scale = TRUE,
    center = TRUE,
    BPPARAM = cores
  ) |>
    aggregate(bin = bin_size)

  p <- ce |>
    as_tibble() |>
    ggplot(aes(x = coord, y = mean, fill = factor(.sample, levels = names(bigwigs)), color = factor(.sample, levels = names(bigwigs)))) + # Ordering the color and fill by the order given in the bigwig input:
    # geom_ribbon(aes(ymin = ci_low, ymax = ci_high), linetype = 2, size = .1, alpha = 0.025) +
    geom_line(size = .5) +
    facet_wrap(features ~ ., ncol = ncols, nrow = nrows) +
    xlab("Distance (bp)") +
    ylab("Mean coverage\n(95% CI)") +
    scale_color_nejm() +
    theme(legend.position = "top") +
    guides(
      color = guide_legend(title = "Track"),
      fill = guide_legend(title = "Track")
    )
  return(p)
}

#### This is simple function to calculate a fit for an Euler plot for the eulerr package. Takes a list of GRanges as an input.
plot_euler <- function(gr_list) {
  # Ensure the input is a list of GRanges objects
  if (!all(sapply(gr_list, inherits, "GRanges"))) {
    stop("All elements in the list must be GRanges objects")
  }

  ovl_final_n <- unlist(sapply(seq_along(gr_list), function(i) {
    combs <- combn(gr_list, i, simplify = F)

    combs_ovl_n <- sapply(seq_along(combs), function(i) {
      suppressWarnings(length(purrr::reduce(combs[[i]], join_overlap_inner)))
    })

    names(combs_ovl_n) <- sapply(seq_along(combs), function(i) {
      paste0(names(combs[[i]]), collapse = "&")
    })
    combs_ovl_n
  }))
  fit <- euler(ovl_final_n)
  return(fit)
}
