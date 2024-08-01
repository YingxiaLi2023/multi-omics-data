########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Results")

########################################################

##### load data 
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
load("./rda_files/resultsumsum.RData")


head(resultsumsum)


ranksall <- resultsumsum

ibrierinds <- grep("ibrier_", names(resultsumsum))
datasets <- unique(resultsumsum$dat)
for(i in seq(along=ibrierinds)) {
  for(j in seq(along=datasets)) {
    ranksall[which(ranksall$dat==datasets[j]),ibrierinds[i]] <- rank(ranksall[which(ranksall$dat==datasets[j]),ibrierinds[i]])
  }
}

cindexinds <- grep("cindex_", names(resultsumsum))
datasets <- unique(resultsumsum$dat)
for(i in seq(along=cindexinds)) {
  for(j in seq(along=datasets)) {
    ranksall[which(ranksall$dat==datasets[j]),cindexinds[i]] <- rank(-ranksall[which(ranksall$dat==datasets[j]),cindexinds[i]])
  }
}


ranksall <- ranksall[ranksall$comb %in% c("rna", "mirna", "methy", "mutation", "cnv", "rna_mirna_methy_mutation_cnv"),]


library(tidyr)             # Load just the tidyr package

# Now you can use pivot_longer directly
ranksall_long <- ranksall %>%
  pivot_longer(
    cols = starts_with("ibrier_") | starts_with("cindex_"),
    names_to = c("metric", "method"),
    names_pattern = "(.*)_(.*)",
    values_to = "rank"
  )


ranksall_long <- as.data.frame(ranksall_long)

ranksall_long$comb <- factor(ranksall_long$comb, levels=c("rna", "mirna", "cnv", "methy", "mutation", "rna_mirna_methy_mutation_cnv"))
levels(ranksall_long$comb) <- c("rna", "mirna", "cnv", "met", "mut", "all (rna, mirna, cnv, met, mut)")

ranksall_long$method <- factor(ranksall_long$method, levels=c("rf", "bf", "lasso", "ipflasso", "prioritylasso"))
levels(ranksall_long$method)[1] <- "rsf"



ns <- c(382, 735, 191, 106, 443, 419, 159, 426, 418, 124, 126, 249, 295, 405)

n_es <- c(103, 72, 17, 37, 152, 77, 35, 101, 132, 52, 38, 62, 62, 38)


ranksall_long$n <- ranksall_long$n_e <- NA

for(i in seq(along=datasets)) {
  ranksall_long$n[ranksall_long$dat==datasets[i]] <- ns[i]
  ranksall_long$n_e[ranksall_long$dat==datasets[i]] <- n_es[i]
}




widths <- 10; heights <- 13

theme_common <- theme(axis.title=element_text(size = 15),
                      strip.text = element_text(size = 14),
                      axis.text = element_text(size = 12, color="black"))




p <- ggplot(ranksall_long[ranksall_long$metric=="ibrier" & ranksall_long$comb %in% c("rna", "mirna", "cnv"),], aes(x=n_e, y=rank)) +
  geom_point() + facet_wrap(~method+comb, ncol=3) + theme_bw() +
  geom_smooth(method = "loess", se = FALSE) + ylab("Ranks") +
  theme_common
p

ggsave(file="./figures/figureS8.png", width=widths, height=heights)



p <- ggplot(ranksall_long[ranksall_long$metric=="ibrier" & ranksall_long$comb %in% c("met", "mut", "all (rna, mirna, cnv, met, mut)"),], aes(x=n_e, y=rank)) +
  geom_point() + facet_wrap(~method+comb, ncol=3) + theme_bw() +
  geom_smooth(method = "loess", se = FALSE) + ylab("Ranks") +
  theme_common
p

ggsave(file="./figures/figureS9.png", width=widths, height=heights)





p <- ggplot(ranksall_long[ranksall_long$metric=="cindex" & ranksall_long$comb %in% c("rna", "mirna", "cnv"),], aes(x=n_e, y=rank)) +
  geom_point() + facet_wrap(~method+comb, ncol=3) + theme_bw() +
  geom_smooth(method = "loess", se = FALSE) + ylab("Ranks") +
  theme_common
p

ggsave(file="./figures/figureS10.png", width=widths, height=heights)



p <- ggplot(ranksall_long[ranksall_long$metric=="cindex" & ranksall_long$comb %in% c("met", "mut", "all (rna, mirna, cnv, met, mut)"),], aes(x=n_e, y=rank)) +
  geom_point() + facet_wrap(~method+comb, ncol=3) + theme_bw() +
  geom_smooth(method = "loess", se = FALSE) + ylab("Ranks") +
  theme_common
p

ggsave(file="./figures/figureS11.png", width=widths, height=heights)








# Split data into subgroups based on combinations of 'comb', 'metric', and 'method'
split_data <- split(ranksall_long, list(ranksall_long$comb, ranksall_long$metric, ranksall_long$method))

# Function to perform linear regression and extract coefficient and p-value
regression_analysis <- function(data) {
  model <- lm(rank ~ n_e, data = data)
  summary_model <- summary(model)
  coef_info <- summary_model$coefficients
  # Return coefficient and p-value for 'n_e' (second row of the coefficients matrix)
  if (nrow(coef_info) > 1) {  # Check if 'n_e' is present in the model
    return(c(coef = coef_info["n_e", "Estimate"], pval = coef_info["n_e", "Pr(>|t|)"]))
  } else {
    return(c(coef = NA, pval = NA))
  }
}

# Apply regression analysis to each subgroup
reg_results <- lapply(split_data, regression_analysis)

# Convert results to a data frame
results_df <- do.call(rbind, reg_results)
row.names(results_df) <- NULL  # Clean up row names for merging later

# Add grouping identifiers back to the results data frame
group_info <- do.call(rbind, lapply(names(reg_results), function(x) strsplit(x, "\\.")))
group_df <- as.data.frame(matrix(unlist(group_info), ncol = 3, byrow = TRUE),
                          stringsAsFactors = FALSE)
colnames(group_df) <- c("comb", "metric", "method")


result_lm <- cbind(group_df, results_df)

result_lm$pval_adj <- NA
result_lm$pval_adj[result_lm$metric=="ibrier"] <- p.adjust(result_lm$pval[result_lm$metric=="ibrier"], method = "holm")
result_lm$pval_adj[result_lm$metric=="cindex"] <- p.adjust(result_lm$pval[result_lm$metric=="cindex"], method = "holm")


tempobj <- result_lm[result_lm$pval < 0.05,]
tempobj[order(tempobj$comb),]

result_lm[result_lm$pval_adj < 0.05,]


result_lm[result_lm$metric=="ibrier",][order(result_lm[result_lm$metric=="ibrier",]$coef),]

result_lm[result_lm$metric=="cindex",][order(result_lm[result_lm$metric=="cindex",]$coef),]


