library(clustAnalytics)
library(pbapply)

# This script contains the code to perform the experiments of the single-level
# SBM networks. 


ground_truth_sbm <- function(blocks=c(40,25,25,10)){
    indices <- 1:length(blocks)
    unlist(mapply(rep, indices, blocks, SIMPLIFY=FALSE))
}

sbm_const <- function(block.sizes=c(25,25,25,25), p_in=0.1, p_out=0.001){
    #sbm with constant p_in and p_out probabilities of inner and outer edges
    n_blocks <- length(block.sizes)
    pm <- matrix(nrow=n_blocks, ncol=n_blocks, data=p_out)
    diag(pm) <- p_in
    g_sbm <- sample_sbm(sum(block.sizes), pref.matrix=pm, block.sizes=block.sizes)
}



scores_sbm <- function(block.sizes, p_in=0.01, p_out=0.001 , ...){
    g <- sbm_const(p_in=p_in, p_out=p_out, block.sizes=block.sizes) %>%
        set_edge_attr("weight", value=1)
    com <- ground_truth_sbm(block.sizes)
    scoring_functions(g, com, ...)
}

merge_scores <- function(block_sizes_list, p_in=0.2, p_out=0.01 , ...){
    scores_list <- pblapply(block_sizes_list, scores_sbm, p_in=p_in, p_out=p_out, ...)
    do.call(rbind, scores_list)
    do.call(rbind, scores_list)
}


blocks_uniform_distribution <- function(lower_bound=10, upper_bound=50, n=10){
    sample(lower_bound:upper_bound, 10)
}


uniform_dist_blocks <- replicate(10000, blocks_uniform_distribution(), simplify=FALSE)



#####################
# Power law blocks
sample_power_law <- function(x_min=1, x_max=10, beta, n){
    #samples according to power law distribution f(x) = x^(-beta)
    y_min <- x_max^(-beta)
    y_max <- x_min^(-beta)
    s <- runif(n, min=y_min, max=y_max)
    s^(-1/beta)
}

sample_block_sizes <- function(n, prob){
    n_com <- length(prob)
    s <- sample(1:n_com, n, prob=prob, replace=TRUE)
    return(as.vector(table(s)))
}

power_law_blocks <- function(network_size=300, n_blocks=15, beta=1.5){
    probs <- sample_power_law(beta=beta, n=n_blocks) 
    # probabilities are not normalized (i.e. they don't sum 1), but the sampling function takes care of that anyway
    sample_block_sizes(network_size, probs)
}

power_law_blocks_varying_size <- function(n=1000, n_blocks_min=5, n_blocks_max=25){
    f <- function(x) replicate(n, power_law_blocks(n_blocks=x), simplify=FALSE)
    l <- lapply(n_blocks_min:n_blocks_max, f)
    cbind(unlist(l, recursive=FALSE))
}


#################


block_sizes_list <- list(rep(10, times=20),
                         rep(25, times=8),
                         rep(50, times=4),
                         rep(100, times=2),
                         c(rep(10, times=10), rep(25, times=4)),
                         c(rep(10, times=5), rep(25, times=2), rep(50, times=2)))



block_sizes_list <- rep(block_sizes_list, times=100) #we generate some duplicates, to have multiple samples of each block structure
set.seed(1)
power_law_block_sizes <- power_law_blocks_varying_size()
scores_table_local <- merge_scores(power_law_block_sizes, type="local")


scores_table_global <- merge_scores(power_law_block_sizes, type="global")
save(scores_table_local, scores_table_global, file="scores_table.RData")



