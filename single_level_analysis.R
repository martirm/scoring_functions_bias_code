# In this script we load the results computed in "single_level_computation.R" 
# and build correlations table to analyze them. This code generates table 3 in 
# the article.

load("scores_table.RData")

cor_matrix_local <- cor(scores_table_local[is.finite(scores_table_local[,"density ratio"]),], use="pairwise.complete.obs")
cor_matrix_global <- cor(scores_table_global[is.finite(scores_table_global[,"density ratio"]),], use="pairwise.complete.obs")

# These are the local and global scores that we will study
scores_local <- c("internal density","edges inside","FOMD","expansion",
                  "cut ratio","conductance", "norm cut", "max ODF","average ODF",
                  "density ratio") 
scores_global <- c("modularity", "coverage", "global density ratio")

# We now build the global and local parts of the correlation table separately
sel_global <- cor_matrix_global[c("n_clusters", "mean_cluster_size"), c(scores_local,scores_global)]
sel_local <- cor_matrix_local[c("size"), scores_local]

# And join them
combined_table <- rbind(sel_global, sel_local)
rownames(combined_table) <- c("# clusters", "mean cluster size", "cluster size")
combined_table[3,c("modularity", "coverage", "global density ratio")] <- NaN #These are global metrics by definition, so we leave these cells empty


# Now use xtable to get the latex code for the table.
library(xtable)
xtable(t(combined_table), digits=4)


