library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

load("multi_level_data.RData")
load("multi_level_BA_data.RData")

scale_colour_Publication <- function(...){
    library(scales)
    discrete_scale("colour","Publication",manual_pal(values = c("#377eb8","#e41a1c","#4daf4a","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
    
}

plot_multi_level_score <- function(score, data, invert_y=FALSE){
    #invert_y: if TRUE, invert Y axis, but only on scores where lower is better.
    score_str <- score
    score <- ensym(score)
    ggplot(data, aes(x = lambda, y = !!score, color=level)) +
        geom_point() + 
        stat_smooth(se = FALSE, size = 1) +
        xlab(expression(lambda)) +
        theme_classic2() + grids(linetype = "dashed") +
        scale_colour_Publication() +
        {
            if (invert_y & score_str %in% reverse_y_scores)
                scale_y_reverse()
        }

}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, 
                                       position = c("bottom", "right"),
                                       title="aaaaa") {
    #source: https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position) + ggtitle(title))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
}


selected_scores <- c("modularity", "coverage", "density ratio", "conductance", "TPR", "norm cut", "expansion", "internal density")
reverse_y_scores <- c("conductance", "norm cut", "expansion")
plot_list <- lapply(selected_scores, plot_multi_level_score, data=multi_level_scores, invert_y=TRUE)
#plot_list <- lapply(selected_scores, plot_multi_level_score, data=multi_level_scores_BA)
do.call(grid_arrange_shared_legend, c(plot_list, list("ncol"=2, "nrow"=4)))


        