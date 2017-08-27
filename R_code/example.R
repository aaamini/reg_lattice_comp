### Source the base code for computing SEM
source("gaussian_projections_base.R") # make sure this file is in your current directory, or add /path/to/ as needed

### Pick a covariance matrix
Sigma <- rbind(
    c(6, 4, -6, -30),
    c(4, 4, -4, -20),
    c(-6, -4, 7, 39),
    c(-30, -20, 39, 234)
)

### Study a single node
jVar <- 2       # which node?
varSet <- c(1)  # which candidate set?
getParentCoefs(jVar, varSet, Sigma)

### Examples: Build a full SEM
neighbourhoodSets <- list(c(2), NULL, c(2,1), c(2, 1, 3))
buildSEM(neighbourhoodSets, Sigma, threshold = TRUE) # use threshold = TRUE to eliminate round-off error

neighbourhoodSets <- list(c(2,3), c(3), NULL, c(2, 1, 3))
buildSEM(neighbourhoodSets, Sigma, threshold = TRUE)

neighbourhoodSets <- list(c(2), c(3), c(4), c(1))
buildSEM(neighbourhoodSets, Sigma, threshold = TRUE)

### Check to see that projecting X(j) onto X(-j) for all j recovers the nonzero entries of the precision matrix
all_nodes <- 1:nrow(Sigma)
neighbourhoodSets <- list(all_nodes[-1], all_nodes[-2], all_nodes[-3], all_nodes[-4])
m <- buildSEM(neighbourhoodSets, Sigma)$coefs
threshold_matrix(m) - diag(diag(m)) # These two matrices should have the
solve(Sigma)                        #  same nonzero entries
