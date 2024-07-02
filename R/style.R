library(BiocParallel)

p <- MulticoreParam(40, progressbar = TRUE)

theme_set(theme_bw() +
            theme(
              line = element_line(linewidth = unit(0.5, "points")),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              panel.border = element_blank(),
              axis.text = element_text(color = "black", size = 10),
              axis.title = element_text(color = "black", size = 10),
              axis.line = element_line(linewidth = unit(0.5, "pt")),
              text = element_text(color = "black"),
              legend.background = element_rect(color = "black", fill = NA, linewidth = 0.5, linetype = "solid")
            ))

bed_colnames = c("seqnames", "start", "end", "name", "score", "strand")

#theme_konsta_minimal <- function (base_size = 12)
#{
#  theme_classic (base_size = base_size)  %+replace%
#    theme (
#      axis.line = element_line(size = .1, color = "black"),
#      axis.line.x = element_line(size = .5, color = "black"),
#      axis.line.y = element_line(size = .5, color = "black"),
#	  axis.text.x = element_text(size = base_size, color = "black"),
#	  axis.text.y = element_text(size = base_size, color = "black"),
#      panel.grid.major.x = element_blank(),
#      panel.grid.minor.y = element_blank(), 
#      panel.grid.major.y = element_blank(),
#      panel.grid.major =element_blank(),
#      panel.grid.minor.x = element_blank(),
#      panel.background = element_rect(fill = NA, color = NA),
#	  plot.title = element_text(size = base_size, color = "black"),
#      legend.background = element_rect(fill=NA, color=NA),
#      legend.key = element_rect(fill=NA, color=NA),
#      strip.text = element_text(size = base_size-1, colour = "black", face = "bold"),
#      strip.background = element_blank(),
#      panel.spacing = unit(.5, "lines"),
#      plot.tag = element_text(face = "bold", size = 12)
#    )
#}

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