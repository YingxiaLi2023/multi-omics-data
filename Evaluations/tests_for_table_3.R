# Test whether the observation that the best combinations featured only few blocks is statistically
# significant.
############################################################################################

# Determine all possible combinations from the five blocks:

blocks <- c("mut", "met", "cnv", "mirna", "rna")

generate_combinations <- function(blocks, length) {
  comb <- expand.grid(rep(list(blocks), length))
  unique(apply(comb, 1, function(row) paste(sort(unique(row)), collapse = " ")))
}

all_combinations <- lapply(1:length(blocks), function(len) generate_combinations(blocks, len))
all_combinations <- unlist(all_combinations)
all_combinations <- unique(all_combinations)



# The numbers of blocks in each combination:

length_combinations <- sapply(all_combinations, function(x) length(strsplit(x, split=" ")[[1]]))
names(length_combinations) <- NULL



# Calculate the probability that the number of datasets where the best performing 
# combination according to the ibrier or index contains only one or two types  
# is at least equal to the observed count:

# ibrier (observed count: 11):
1 - pbinom(q=11-1, size=14, prob=mean(length_combinations <= 2))

# cindex (observed count: 14):
1 - pbinom(q=14-1, size=14, prob=mean(length_combinations <= 2))





# Test whether mRNA and miRNA are overrepresented in the best combinations.
###########################################################################

library("poisbinom")

# RNA:

# ibrier:
1 - ppoisbinom(q=8-1, pp=c(rep(1/5, 4), rep(2/5, 7), rep(3/5, 2), 4/5))

# cindex:
1 - ppoisbinom(q=9-1, pp=c(rep(1/5, 8), rep(2/5, 5), 3/5))


# miRNA:

# ibrier:
1 - ppoisbinom(q=6-1, pp=c(rep(1/5, 4), rep(2/5, 7), rep(3/5, 2), 4/5))

# cindex:
1 - ppoisbinom(q=6-1, pp=c(rep(1/5, 8), rep(2/5, 5), 3/5))


# Adjust for multiple testing:

5*(1 - ppoisbinom(q=9-1, pp=c(rep(1/5, 8), rep(2/5, 5), 3/5)))
