library(ggsci)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(soGGi)

# Date today in YYMMDD format:
dt <- function() {
  today <- Sys.Date()
  return(format(today, format = "%y%m%d"))
}

# Seed for reproducibility
set.seed(123)

# Canonical chromosomes for filtering convinience
chrs <- c(paste0("chr", 1:22), "chrX", "chrY")

# Setting date format for plot prefixes:
today <- Sys.Date()
today <- format(today, format = "%y%m%d")

# Colnames for MACS2 output narrowpeak files
narrowpeak_colnames <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pval", "qval", "peak")

# Column names for AME motif discovery out files
ame_cnames <- c("rank", "motif_db", "motif_ID", "motif_name", "consensus", "pval", "adj_pval", "eval", "tests", "fasta_max", "pos", "neg", "pwm_min", "tp", "pct_tp", "fp", "pct_fp")

# Annotation directory
anndir <- "/data/kzkarttu/annotations"
# anndir <- "/data/kzkarttu/gradu/annotations"

# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://nov2020.archive.ensembl.org")

gencode_annotation <- read_tsv(file.path(anndir, "genes/gencode.v36.annotation.genebody.bed"), col_names = c("seqnames", "start", "end", "ensembl_gene_id_version", "score", "strand", "hgnc_symbol"))
gencode_annotation_tss <- read_tsv(file.path(anndir, "genes/gencode.v36.annotation.TSS.bed"), col_names = c("seqnames", "start", "end", "ensembl_gene_id_version", "score", "strand", "hgnc_symbol"))

gencode_annotation <- gencode_annotation %>% 
  mutate(ensid = str_remove(ensembl_gene_id_version, "\\.\\d+"))
gencode_annotation_tss <- gencode_annotation_tss %>% 
  mutate(ensid = str_remove(ensembl_gene_id_version, "\\.\\d+"))

# if (!exists("repeatmasker")) {
#   repeatmasker <- read_tsv(file.path("/data/kzkarttu/gradu/annotations/repeatMasker/repeatMasker_sorted.bed"), col_names = c("seqnames", "start", "end", "name", "score", "strand")) %>%
#     filter(seqnames %in% chrs) %>% # Only canonical chromosomes (chr1-22,XY)
#     mutate(
#       class = str_extract(name, "[[:alnum:]:punct:]+$"), # Making columns for TE classes and subfamilies
#       subf = str_remove(name, "chr\\d+\\|\\d+\\|\\d+\\||chr[XY]\\|\\d+\\|\\d+\\|"),
#       te_name = str_extract(subf, "[^\\|]+")
#     )
#   repeatmasker_gr <- GRanges(repeatmasker)
# }

# Repeatmasker annotations:
if (!exists("rmsk")) {
  rmsk <- read_tsv(paste0(anndir, "/repeats/repeatMasker_filtered.bed.gz"), col_types = "ciiccc", col_names = c("seqnames", "start", "end", "name", "score", "strand")) %>% # Specifying col_types to prevent scientific notation for start-end coordinates
    filter(seqnames %in% chrs) %>% # Only canonical chromosomes (chr1-22,XY)
    mutate(
      class = str_extract(name, "LINE$|SINE$|LTR$|DNA$|Retroposon$"), # Making columns for TE classes and subfamilies
      family = str_remove_all(name, "^[^\\|]*\\||\\|[:alnum:]+$"),
      subfamily = str_extract(name, "[^\\|]*")
    )
  rmsk_gr <- GRanges(rmsk)
}

# Some dplyr options
options(
  readr.num_columns = 0,
  dplyr.summarise.inform = F
)

# Bedtools path for bedtoolsr package
options(bedtools.path = "/data/miniconda2/envs/gradu_3.8/bin/")

# Fonts and color palettes
# font_add("Helvetica",
#   regular = "/data/kzkarttu/gradu/scripts/R_scripts/gradu/font/Helvetica.ttf", bold = "/data/kzkarttu/gradu/scripts/R_scripts/gradu/font/Helvetica-Bold.ttf",
#   italic = "/data/kzkarttu/gradu/scripts/R_scripts/gradu/font/Helvetica-Oblique.ttf", bolditalic = "/data/kzkarttu/gradu/scripts/R_scripts/gradu/font/Helvetica-BoldOblique.ttf", symbol = NULL
# )

# hg19 to hg38 chain 
chain <- rtracklayer::import.chain("/data/kzkarttu/gradu/annotations/UCSC/hg19ToHg38.over.chain")
chrom_sizes <- "/data/kzkarttu/gradu/annotations/UCSC/hg38.canonical.chrom.sizes"
mask <- "/data/kzkarttu/gradu/annotations/mask/BSgenome_hg38_mask.bed"

# font_add(family = "Arial", regular = "/data/kzkarttu/gradu/scripts/R_scripts/gradu/font/arial.ttf")
# showtext_auto()
# windows()

#theme_set(theme_bw() +
#  theme(
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    strip.background = element_blank(),
#    panel.border = element_rect(size = 0.5, colour = "black"),
#    text = element_text(color = "black", family = "Arial"),
#    axis.text = element_text(color = "black", family = "Arial"),
#    legend.background = element_rect(color = "black", fill = NA, size = 0.5, linetype = "solid")
#  ))

# Color palettes
col_pal <- c(pal_npg(palette = "nrc")(4), "gray40")
col_pal_nejm <- pal_nejm("default")(8)
cluster_pal <- pal_nejm("default")(5)
heat_cols_wes <- wes_palette("Zissou1", n = 5)
heat_pal_wes <- colorRampPalette(c(heat_cols_wes[2], heat_cols_wes[3], heat_cols_wes[5]))(100)
heat_cols_red <- c(scales::muted("blue", 40, 70), "white", scales::muted("red", 40, 100))
heat_cols_orrd <- RColorBrewer::brewer.pal(9, "OrRd")
heat_cols_orteal <- c("#008083", "#259EA1", "#FBAB37", "#F68105", "#FC5802")
#"#036061", 
# color_palette = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

# Function for detaching packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats", "package:graphics", "package:grDevices", "package:utils", "package:datasets", "package:methods", "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0) for (package in package.list) detach(package, character.only = TRUE)
}

render_with_jobs <- function() {
  rstudioapi::verifyAvailable()
  jobs_file <- tempfile(tmpdir = "/tmp", fileext = ".R")
  rmd_to_render <- rstudioapi::selectFile(
    caption = "Choose an Rmd file...",
    filter = "Rmd files (*.Rmd)"
  )
  if (is.null(rmd_to_render)) {
    stop("You must choose an Rmd file to proceed!")
  }
  cat(paste0('rmarkdown::render("', rmd_to_render, '")'), file = jobs_file)
  rstudioapi::jobRunScript(
    path = jobs_file,
    name = basename(rmd_to_render),
    workingDir = getwd()
  )
  # Wait for 5 seconds so that the tmp file isn't deleted *before* the job
  # submission is complete
  Sys.sleep(5)
  file.remove(jobs_file)
}

###########################################################################################################################################
##################################################### Metaplotting function ###############################################################
###########################################################################################################################################

plot_profile <- function(bws, regs, size, rollmean_window = 75, cores = 40) {
  if (is.list(regs) == F) { # Checking if there is only one table of regions or if the regs is a list; if so, run the function without a loop over the list

    regs_cov <- bind_rows(lapply(seq_along(bws), function(i) {
      reg <- regionPlot(bws[i], testRanges = regs, format = "bigwig", distanceAround = reg_size, normalize = T)

      colMeans(assay(reg)) %>%
        as_tibble() %>%
        mutate(
          coord = row_number(),
          te_name = unique(regs$te_name),
          sample = names(bws)[i]
        )
    }))
  } else { # If regs is a list, loop over the regs in parallel
    
    regs_cov <- bind_rows(mclapply(seq_along(regs), mc.cores = cores, function(y) {
      bind_rows(lapply(seq_along(bws), function(i) {
        reg <- regionPlot(bws[i], testRanges = regs[[y]], format = "bigwig", distanceAround = reg_size, normalize = T)

        if (!is.null(regs[[y]]$cluster)) {
          colMeans(assay(reg)) %>%
            as_tibble() %>%
            mutate(
              coord = row_number(),
              te_name = unique(regs[[y]]$te_name),
              sample = names(bws)[i],
              cluster = unique(regs[[y]]$cluster)
            )
        } else {
          colMeans(assay(reg)) %>%
            as_tibble() %>%
            mutate(
              coord = row_number(),
              te_name = unique(regs[[y]]$te_name),
              sample = names(bws)[i]
            )
        }
      }))
    }))
  }

  if (!is.null(regs_cov$cluster)) {
    regs_cov %>%
      group_by(te_name, sample, cluster) %>%
      mutate(sample = factor(sample, levels = unique(names(bws)))) %>%
      mutate(rollmean = zoo::rollmean(log2(value + 1), rollmean_window, fill = NA))
  } else {
    regs_cov %>%
      group_by(te_name, sample) %>%
      mutate(sample = factor(sample, levels = unique(names(bws)))) %>%
      mutate(rollmean = zoo::rollmean(log2(value + 1), rollmean_window, fill = NA))
  }
}

###########################################################################################################################################
############################################## Metaplotting function, with no grouping for TEs ############################################
###########################################################################################################################################

plot_profile_2 <- function(bws, regs, size, rollmean_window = 75, cores = 40) {
  if (is.list(regs) == F) { # Checking if there is only one table of regions or if the regs is a list; if so, run the function without a loop over the list
    
    regs_cov <- bind_rows(lapply(seq_along(bws), function(i) {
      reg <- regionPlot(bws[i], testRanges = regs, format = "bigwig", distanceAround = reg_size, normalize = T)
      
      colMeans(assay(reg)) %>%
        as_tibble() %>%
        mutate(
          coord = row_number(),
          sample = names(bws)[i]
        )
    }))
  } else { # If regs is a list, loop over the regs in parallel
    
    regs_cov <- bind_rows(mclapply(seq_along(regs), mc.cores = cores, function(y) {
      bind_rows(lapply(seq_along(bws), function(i) {
        reg <- regionPlot(bws[i], testRanges = regs[[y]], format = "bigwig", distanceAround = reg_size, normalize = T)
        
        if (!is.null(regs[[y]]$cluster)) {
          colMeans(assay(reg)) %>%
            as_tibble() %>%
            mutate(
              coord = row_number(),
              sample = names(bws)[i],
              cluster = unique(regs[[y]]$cluster)
            )
        } else {
          colMeans(assay(reg)) %>%
            as_tibble() %>%
            mutate(
              coord = row_number(),
              sample = names(bws)[i]
            )
        }
      }))
    }))
  }
  
  if (!is.null(regs_cov$cluster)) {
    regs_cov %>%
      group_by(sample, cluster) %>%
      mutate(sample = factor(sample, levels = unique(names(bws)))) %>%
      mutate(rollmean = zoo::rollmean(value, rollmean_window, fill = NA))
  } else {
    regs_cov %>%
      group_by(sample) %>%
      mutate(sample = factor(sample, levels = unique(names(bws)))) %>%
      mutate(rollmean = zoo::rollmean(value, rollmean_window, fill = NA))
  }
}

####### INPUTS FOR THE FUNCTION ###########
# bams: A named vector of .bam files
# regs: A list of ranges in a data frame/tibble,
# The ranges need to have columns for chromosomes (chr), start position (start), end position (end), name of the region (geneid) and width of the region (width)
# Paired: A boolean vector of the same length as bams, determines if bam is paired or not for the featureCounts function 

run_featurecounts <- function(regs, bams, paired, threads = 60) {
    
  bind_cols(lapply(seq_along(bams), function(i) {
      fc <- featureCounts(files = bams[i], annot.ext = regs, isPairedEnd = paired[i], nthreads = threads, verbose = F, useMetaFeatures = F)
      
      #### FPKM normalization ####
      read_n <- sum(fc$stat[, 2]) / 10e6 # Read count in whole library divided by the scaling factor i.e. 1 000 000
      
      fc <- fc$counts / read_n # RPM normalization
      fc <- fc / (regs$width / 1000) # RPKM
      
      fc <- as_tibble(fc)
      colnames(fc) <- names(bams)[i]
      fc
    }))
}







