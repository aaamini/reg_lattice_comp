source("lattice_comp.R")
library(igraph)
library(matrixStats)

dag_cov <- function(B,Omega){
  d <- dim(B)[1]
  R <- eye(d) - B
  Gamma <- drop0(zapsmall(R %*% solve(Omega) %*% t(R)))
  Sig <- drop0(zapsmall(solve(Gamma)))
  return( list(Gamma=Gamma, Sig=Sig))
}

data_file <- 'Bet_example_3a.RData'

if (file.exists(data_file)) {
  load(data_file)
  d <- dim(Bet)[1]
 
  out <- dag_cov(Bet,Omega)
  Gamma <- out$Gamma
  Sig <- out$Sig
  
} else {
  d <- 10
  Bet <- rsparsematrix(d,d,density=.3, symmetric = T)
  Bet[lower.tri(Bet,diag=T)] <- 0
  Omega <- eye(d)
  
  out <- dag_cov(Bet,Omega)
  Gamma <- out$Gamma
  Sig <- out$Sig
  
  all_nodes <- 1:d
  pcg <- graph_from_adjacency_matrix(as(Gamma,"nMatrix")*1, mode = "undirected", diag = F)
  layout <- layout.fruchterman.reingold(pcg)
  V(pcg)$color <- "#FFFF80"
  V(pcg)$label <- 1:d
  V(pcg)$x <- layout[,1]
  V(pcg)$y <- layout[,2]
  
  dag <- graph_from_adjacency_matrix(as(Bet,"nMatrix")*1, mode = "directed", diag = F)
  V(dag)$x <- V(pcg)$x
  V(dag)$y <- V(pcg)$y
  V(dag)$color <- "#FFFF80"
  V(dag)$label <- 1:d
  E(dag)$arrow.size <- .5
  
  #save(Bet, Omega, pcg, dag, file="Bet_example_3a.RData")
}


A <- get.adjacency(pcg)
image(A)

A[upper.tri(A)] <- 0
A[as.matrix((abs(A) > 0) & lower.tri(A))] <- rnorm(nnzero(A)/2)
A <- drop0(A)
A <- A + t(A)

non_psd <- TRUE
while (non_psd){
  eps <- runif(1)
  new_G <- diag(1,d) + eps*A
  eig_vals <- eigen(new_G)$values
  non_psd <- min(eig_vals) < 0
  cat('.')
}
  
Gamma <- new_G
Sig <- solve(Gamma)
