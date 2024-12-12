# Rmarkdown options
knitr::opts_chunk$set(
  cache = FALSE, warning = FALSE,
  message = FALSE, cache.lazy = FALSE,
  fig.height = 4, fig.width = 6,
  fig.asp = 0.618, fig.show = "show",
  out.width = "70%", fig.align = "center"
)

# Some dplyr options
options(
  readr.num_columns = 0,
  dplyr.summarise.inform = F
)

# GGplot theme, adapted from cowplot package:
theme_cow <- function(font_size = 14, font_family = "", line_size = 0.5,
                      rel_small = 12 / 14, rel_tiny = 11 / 14, rel_large = 16 / 14) {
  half_line <- font_size / 2
  small_size <- rel_small * font_size
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      line = element_line(
        color = "black", linewidth = line_size,
        linetype = 1, lineend = "butt"
      ), rect = element_rect(
        fill = NA,
        color = NA, linewidth = line_size, linetype = 1
      ),
      text = element_text(
        family = font_family, face = "plain",
        color = "black", size = font_size, hjust = 0.5,
        vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(),
        debug = FALSE
      ), axis.line = element_line(
        color = "black",
        linewidth = line_size, lineend = "square"
      ),
      axis.line.x = NULL, axis.line.y = NULL, axis.text = element_text(
        color = "black",
        size = small_size
      ), axis.text.x = element_text(
        margin = margin(t = small_size / 4),
        vjust = 1
      ), axis.text.x.top = element_text(
        margin = margin(b = small_size / 4),
        vjust = 0
      ), axis.text.y = element_text(
        margin = margin(r = small_size / 4),
        hjust = 1
      ), axis.text.y.right = element_text(
        margin = margin(l = small_size / 4),
        hjust = 0
      ), axis.ticks = element_line(
        color = "black",
        linewidth = line_size
      ), axis.ticks.length = unit(
        half_line / 2,
        "pt"
      ), axis.title.x = element_text(
        margin = margin(t = half_line / 2),
        vjust = 1
      ), axis.title.x.top = element_text(
        margin = margin(b = half_line / 2),
        vjust = 0
      ), axis.title.y = element_text(
        angle = 90,
        margin = margin(r = half_line / 2), vjust = 1
      ),
      axis.title.y.right = element_text(
        angle = -90, margin = margin(l = half_line / 2),
        vjust = 0
      ), legend.background = element_blank(),
      legend.spacing = unit(font_size, "pt"), legend.spacing.x = NULL,
      legend.spacing.y = NULL,
      legend.margin = margin(0, 0, 0, 0),
      legend.key = element_blank(),
      legend.key.size = unit(1.1 * font_size, "pt"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(rel_small)),
      legend.text.align = NULL,
      legend.title = element_text(hjust = 0),
      legend.title.align = NULL,
      legend.position = "right",
      legend.direction = NULL,
      legend.justification = "center", legend.box = NULL, legend.box.margin = margin(0, 0, 0, 0),
      legend.box.background = element_blank(),
      legend.box.spacing = unit(font_size, "pt"),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = NULL, panel.grid.minor = NULL,
      panel.grid.major.x = NULL, panel.grid.major.y = NULL,
      panel.grid.minor.x = NULL, panel.grid.minor.y = NULL,
      panel.spacing = unit(half_line, "pt"), panel.spacing.x = NULL,
      panel.spacing.y = NULL, panel.ontop = FALSE, strip.background = element_rect(fill = "white"),
      strip.text = element_text(
        size = rel(rel_small),
        margin = margin(
          half_line / 2, half_line / 2, half_line / 2,
          half_line / 2
        )
      ), strip.text.x = NULL, strip.text.y = element_text(angle = -90),
      strip.placement = "inside", strip.placement.x = NULL,
      strip.placement.y = NULL, strip.switch.pad.grid = unit(
        half_line / 2,
        "pt"
      ), strip.switch.pad.wrap = unit(
        half_line / 2,
        "pt"
      ), plot.background = element_blank(),
      plot.title = element_text(
        face = "bold",
        size = rel(rel_large), hjust = .5, vjust = 1,
        margin = margin(b = half_line),
      ), plot.subtitle = element_text(
        size = rel(rel_small),
        hjust = 0, vjust = 1, margin = margin(b = half_line)
      ),
      plot.caption = element_text(
        size = rel(rel_tiny),
        hjust = 1, vjust = 1, margin = margin(t = half_line)
      ),
      plot.tag = element_text(
        face = "bold", hjust = 0,
        vjust = 0.7
      ), plot.tag.position = c(0, 1), plot.margin = margin(
        half_line,
        half_line, half_line, half_line
      ), complete = TRUE
    )
}

ggplot2::theme_set(theme_cow())

# theme_set(theme_bw() +
#   theme(
#     line = element_line(linewidth = unit(0.5, "points")),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank(),
#     axis.text = element_text(color = "black", size = 10),
#     axis.title = element_text(color = "black", size = 10),
#     axis.line = element_line(linewidth = unit(0.5, "pt")),
#     text = element_text(color = "black"),
#     # legend.background = element_rect(color = "black", fill = NA, linewidth = 0.5, linetype = "solid")
#     legend.background = element_blank(),
#     legend.position = "bottom",
#     legend.justification = "center",
#     legend.key.spacing = unit(5.5, "points"),
#     legend.box.spacing = unit(11, "points"),
#     legend.ticks.length = " "
#   ))

# New default theme, inspired by ArchR
# ggplot2::theme_set(theme_bw() +
#            theme(
#              text = element_text(family = "sans"),
#              axis.text = element_text(color = "black", size = 10),
#              axis.title = element_text(color = "black", size = 10),
#              title = element_text(color = "black", size = 10),
#              plot.margin = unit(c(1, 1, 1, 1), "cm"),
#              panel.background = element_rect(fill = "transparent", color = NA),
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              axis.line = element_line(linewidth = unit(0.5, "pt"), size = (4 / 3) * .5 * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
#              panel.border = element_blank(),
#              # panel.border = element_rect(fill = NA, color = "black", size = (4 / 3) * .5 * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
#              axis.ticks.length = unit(0.1, "cm"),
#              axis.ticks = element_line(color = "black", size = .5 * (4 / 3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
#              legend.key = element_rect(fill = "transparent", colour = NA),
#              legend.text = element_text(color = "black", size = 8),
#              legend.box.background = element_rect(color = NA),
#              legend.position = "bottom",
#              strip.text = element_text(size = 10, color = "black")
#            ))


# theme_konsta_minimal <- function (base_size = 12)
# {
#  theme_classic (base_size = base_size)  %+replace%
#    theme (
#      axis.line = element_line(size = .1, color = "black"),
#      axis.line.x = element_line(size = .5, color = "black"),
#      axis.line.y = element_line(size = .5, color = "black"),
# 	  axis.text.x = element_text(size = base_size, color = "black"),
# 	  axis.text.y = element_text(size = base_size, color = "black"),
#      panel.grid.major.x = element_blank(),
#      panel.grid.minor.y = element_blank(),
#      panel.grid.major.y = element_blank(),
#      panel.grid.major =element_blank(),
#      panel.grid.minor.x = element_blank(),
#      panel.background = element_rect(fill = NA, color = NA),
# 	  plot.title = element_text(size = base_size, color = "black"),
#      legend.background = element_rect(fill=NA, color=NA),
#      legend.key = element_rect(fill=NA, color=NA),
#      strip.text = element_text(size = base_size-1, colour = "black", face = "bold"),
#      strip.background = element_blank(),
#      panel.spacing = unit(.5, "lines"),
#      plot.tag = element_text(face = "bold", size = 12)
#    )
# }

# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#
#   numPlots = length(plots)
#
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#
#   if (numPlots==1) {
#     print(plots[[1]])
#
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }
