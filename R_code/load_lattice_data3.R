source("lattice_comp.R")
library(igraph)
library(matrixStats)

random_invcov <- function(A) {
  #A[upper.tri(A)] <- 0
  #A[as.matrix((abs(A) > 0) & lower.tri(A))] <- rnorm(nnzero(A))
  #A <- drop0(A)
  #A <- A + t(A)
  
  non_psd <- TRUE
  while (non_psd) {
    eps <- runif(1)
    Gamma <- diag(1,d) + eps*A
    eig_vals <- eigen(Gamma)$values
    non_psd <- min(eig_vals) < 0
    cat('.')
  }
  list(Gamma=Gamma, eps=eps) 
}


data_file <- 'Sig_example_4a.RData'
#data_file <- 'Sig_example_4.RData'


if (file.exists(data_file)) {
  load(data_file)
  d <- dim(Sig)[1]

} else {
  d <- 15
  A <- rsparsematrix(d,d,density=.15, symmetric = T)
  diag(A) <- 0
  
  
  Gamma <- random_invcov(A)$Gamma
  Sig <- solve(Gamma)
  
  all_nodes <- 1:d
  pcg <- graph_from_adjacency_matrix(as(Gamma,"nMatrix")*1, mode = "undirected", diag = F)
  layout <- layout.fruchterman.reingold(pcg)
  #layout <- layout.drl(pcg)
  V(pcg)$color <- "#FFFF80"
  V(pcg)$label <- 1:d
  V(pcg)$x <- layout[,1]
  V(pcg)$y <- layout[,2]
  
  V(pcg)$size <- 20
  V(pcg)$label.cex <- 0.9
  plot(pcg)
  save(A, Gamma, Sig, d, pcg, all_nodes, file=data_file)
}



