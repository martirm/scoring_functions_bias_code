library(clustAnalytics)
library(pbapply)

# This script contains the code to perform the experiments of the multi-level
# preferential attachment graph based on the Barabási-Albert model. The output 
# is saved at "multi_level_BA_data.RData"
# to be plotted in "multi_level_analysis.R".

build_affinity_matrix <- function(block_sizes_UL = c(2,2,2,2), p1, p2, p3){
    n <- sum(block_sizes_UL)
    M <- matrix(nrow=n, ncol=n, p3)
    i <- 1
    for (j in block_sizes_UL){
        rg <- i:(i+j-1)
        M[rg,rg] <- p2
        i <- i+j
    }
    diag(M) <- p1
    return(M)
}


BA_multi_level <- function(block_sizes_UL=c(2,2,2,2), label_probs_LL, 
                           b1, b2, b3, m, t_max){
    # generates two-level sbm of the given block sizes and probabilites
    # block_sizes1: higher level (number of lower level communities it includes, not vertices)
    # block_sizes2: lower_level
    # b1: probability of inner edges in the lower level of communities
    # b2: inner edges on the higher level, but outer edges on the lower
    # b3: outer edges on both
    
    B <- build_affinity_matrix(block_sizes_UL, b1, b2, b3)
    p <- label_probs_LL
    #n <- sum(block_sizes_UL*block_sizes_LL)
    #m <- expected_number_edges(B, block_sizes_LL)
    #edges_per_step <- m/n
    barabasi_albert_blocks(m=m,p=p, B=B, t_max=t_max, 
                           sample_with_replacement = TRUE, type="block_first")
    
    
}

labels_upper_level <- function(labels_LL, block_sizes_UL){
    labels_UL <- labels_LL
    x <- 1:sum(block_sizes_UL) #x[l] is the upper level membership of lower level label l
    j <- 1
    for (i in 1:length(block_sizes_UL)){
        x[j:(j+block_sizes_UL[i]-1)] <- i
        j <- j + block_sizes_UL[i]
    }
    return(x[labels_LL])
}

BA_multi_level_scores_table <- function(block_sizes_UL=c(2,2,2,2), label_probs_LL=label_probs_LL, 
                                         b1, b3, m, t_max=200, length.out=100){
    #p2 will vary, p1 and p3 are fixed
    # length.out: the number of values p2 will take (which will determine the "resolution" of the 
    # plot)
    

    lambda_seq <- seq(from=0, to=1, length.out=length.out)
    sample_graph <- function(lambda){
        g <- BA_multi_level(block_sizes_UL, label_probs_LL, b1, b3+lambda*(b1-b3), b3, m, t_max)
        scores_lower <- scoring_functions(g, V(g)$label, type="global", weighted=FALSE)
        scores_upper <- scoring_functions(g, labels_upper_level(V(g)$label, block_sizes_UL), type="global", weighted=FALSE)
        list("lower"=scores_lower, "upper"=scores_upper)
    }
    all_scores <- lapply(lambda_seq, sample_graph)
    scores_lower <- cbind(as.data.frame(do.call(rbind, lapply(all_scores, `[[`, 1))), "lambda"=lambda_seq, "level"="lower")
    scores_upper <- cbind(as.data.frame(do.call(rbind, lapply(all_scores, `[[`, 2))), "lambda"=lambda_seq, "level"="upper")
    
    rbind(scores_lower, scores_upper)
}

set.seed(1)
multi_level_scores_BA <- BA_multi_level_scores_table(block_sizes_UL=c(2,2,2,2), label_probs_LL=rep(1/8, 8), 
                                                   b1=0.2, b3=.01, m=4, t_max=300, length.out=100)
save(multi_level_scores_BA, file="multi_level_BA_data.RData")



