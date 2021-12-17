load("scores_table.RData")


cor_matrix_local <- cor(scores_table_local[is.finite(scores_table_local[,"density ratio"]),], use="pairwise.complete.obs")
cor_matrix_global <- cor(scores_table_global[is.finite(scores_table_global[,"density ratio"]),], use="pairwise.complete.obs")

#Add back TPR once it's implemented in the rewire package!!!
scores_local <- c("internal density","edges inside","FOMD","expansion",
                  "cut ratio","conductance", "norm cut", "max ODF","average ODF",
                  "density ratio") 
scores_global <- c("modularity", "coverage", "global density ratio")

sel_global <- cor_matrix_global[c("n_clusters", "mean_cluster_size"), c(scores_local,scores_global)]
sel_local <- cor_matrix_local[c("size"), scores_local]
combined_table <- rbind(sel_global, sel_local)
rownames(combined_table) <- c("# clusters", "mean cluster size", "cluster size")
combined_table[3,c("modularity", "coverage", "global density ratio")] <- NaN

library(xtable)
xtable(t(combined_table), digits=4)

##############
#linearMod <- lm(size ~ TPR, data=scores_table_local)
