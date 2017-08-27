### Auxiliary function
ones <- function(m,n) matrix(rep(1,m*n),nrow=m)
zeros <- function(m,n) matrix(rep(0,m*n),nrow=m)
eye <- function(n) diag(rep(1,n))

DEFAULT.TOL <- 1e-9

#### Lattice related 
library(Matrix)

threshold <- function(m, tol){
  m[abs(m) < tol] <- 0
  
  m
}

getParentCoefs2 <- function(j, S, sigma, threshold=F, tol=DEFAULT.TOL){
  d <- dim(sigma)[1]
  beta <- rep(0,d)
  beta <- as(t(beta),"sparseMatrix")
  colnames(beta) <- 1:d
  
  if( !is.null(S) && (length(S) > 0) ){
    sigmaXY <- sigma[j, S]
    sigmaYY <- sigma[S, S]
    coefs <- sigmaXY %*% solve(sigmaYY)
  
    # if (threshold) coefs <- zapsmall(coefs)
    if (threshold) coefs <- threshold(coefs, tol)
    beta[S] <- coefs
  }
  
  list(coefs = beta, parents = S, node = j)
}

computeLattice <- function(j, S, Sig, all_nodes, restricted_set=NULL, should_sort=F) {
  supp <- function(c) all_nodes[as.vector(as(c,"nMatrix"))] 
  
  if (is.null(restricted_set))  restricted_set <- all_nodes
  aux <- setdiff(restricted_set, c(j,S))
  M <- S
  c <- getParentCoefs2(j, S, Sig, threshold=T)$coefs
  m <- supp(c)
  while ( length(aux) ) {
    k <- aux[1]
    propose <- union(k,M)
    aux <- setdiff(aux,k)
    c <- getParentCoefs2(j, propose, Sig, threshold=T)$coefs
    if ( setequal(supp(c), m) ) {
      M <- union(M,k)
    }
  }
  
  
  sort_or_not <- function(x) if (should_sort) sort(x) else x
  m=sort_or_not(m)
  M=sort_or_not(M)
  S=sort_or_not(S)
  
  class(m) <- 'set'
  class(M) <- 'set'
  class(S) <- 'set'
  
  lattice <- list(m=m, M=M,num_elms = 2^length(setdiff(M,m)), S=S)
  class(lattice) <- "nbhdLattice"
  return ( lattice )
}

print.set <- function(set, width=NA){
  set_string_repr <- paste(set,collapse=",")
  
  if (!is.na(width)) {
    n <- nchar(set_string_repr)+2
    if (n < width) {
      cat( paste(replicate(width-n," "),collapse="") )
    }
  }
  cat('{')
  cat(set_string_repr)
  cat('}')
}
print.nbhdLattice <- function(lattice, width=NA, VERB=1){
  
  print(lattice$m, width=width)
  cat(' -- ')
  print(lattice$M, width=width)
  #cat(sprintf('   (%d sets)\n',lattice$num_elms))
  temp_str <- sprintf('(%g sets)',lattice$num_elms)
  cat(sprintf('%10s',temp_str)) 
  if (VERB > 1) {
    cat(sprintf('%7s','S = '))
    print(structure(lattice$S,class="set"))  # ugly hack!  have to implement as(...,"set")
  }
}

genLattice <- function(lattice) {
  setdiff(lattice$M, lattice$m)
}


trim_collection_A_by_B <- function(A, B) {
  orig_len <- length(A)
  for (k in 1:orig_len) {
    if ( is_set_in_collection(A[[k]], B) )  A[[k]] <- NA
  }
  A <- A[!is.na(A)]
  new_len <- length(A)
  if (orig_len != new_len)  cat(sprintf('trimmed A: %d -> %d \n',orig_len, new_len))
  
  return(A)
}
#A <- list(c(1,2,3),c(1,4),c(5,6))
#B <- list(c(1,4),c(2,3))
#trim_collection_A_by_B(A,B)
trim_collection_A_by_func <- function(A, func, VERB=1) {
  orig_len <- length(A)
  for (k in 1:orig_len) {
    if ( func(A[[k]]) )  A[[k]] <- NA
  }
  A <- A[!is.na(A)]
  new_len <- length(A)
  if (orig_len != new_len && VERB > 2)  cat(sprintf('trimmed A: %d -> %d \n',orig_len, new_len))
  
  return(A)
}


is_set_in_collection <- function(set, collection) {
  any(sapply(collection, function(x) setequal(set,x)))
}

is_set_in_lattices <- function(set, lattices) {
  any(sapply(lattices, function(lat) isInLattice(set,lat)))
}

is.subset <- function(set,superset) {
  return (length(setdiff(set, superset)) == 0)
}
isInLattice <- function(set,lattice){
  return (is.subset(lattice$m, set) && is.subset(set, lattice$M))
}

computeAllSubordinateLattices <- function(j, Sig, all_nodes, root_S = NA, print.width=NA, VERB=1, max.iter=1e4) { 
  if (any(is.na(root_S))) {
    root_S = all_nodes[-j]
  }
  S_candidates <- list()
  S_candidates[[1]] <- root_S 
  #m_idx <- 1
  #minimals <- list() # This is not needed
  
  lattices <- list()
  lat_idx <- 1
  total_num_of_sets <- 0
  for (t in 1:max.iter) {
    
    S <- S_candidates[[1]]
    S_candidates[[1]] <- NULL
   
    lattice <- computeLattice(j, S, Sig, all_nodes)
    print(lattice, width=print.width, VERB=VERB)
    
    lattices[[lat_idx]] <- lattice
    lat_idx <- lat_idx + 1
    total_num_of_sets <- total_num_of_sets + lattice$num_elms
    
    m <- lattice$m
    m_len <- length(m)
    #minimals[[m_idx]] <- m
    #m_idx <- m_idx + 1
    
    if (m_len > 0) {
      potentials <- combn(m, m_len-1, simplify = F)
      
      # S_candidates <- trim_collection_A_by_B( unique( append(S_candidates, potentials) ), minimals )
      S_candidates <- unique( append(S_candidates, potentials) )
      #trim_func <- function(set) {  is_set_in_collection(set, minimals) || is_set_in_lattices(set,lattices) }
      trim_func <- function(set) {  is_set_in_lattices(set,lattices) }
      S_candidates <- trim_collection_A_by_func( S_candidates, trim_func, VERB=VERB)
      
    }
    # if (VERB > 1) {
    #   cat(sprintf('%7s','S = '))
    #   print(structure(S,class="set"))  # ugly hack!  have to complement as(...,"set")
    # }
    cat('\n')
    if (length(S_candidates) == 0) break
  }
  
  
  cat(sprintf("\nTotal number of sets covered = %d\n",total_num_of_sets))
  return(lattices)
}

#verifyCondOrth <- function(A, S, B, Sig, all_nodes, threshold=T, tol=DEFAULT.TOL) {
A_and_B_proj_seperated_by_S <- function(A, S, B, Sig, threshold=T, tol=DEFAULT.TOL) {
  # Test whether P_A P_S^\perp P_B = 0
  if ( !is.null(S) && (length(S) > 0) )
    temp <- Sig[B,A] - Sig[B,S] %*% solve(Sig[S,S]) %*% Sig[S,A]
  else
    temp <- Sig[B,A] 
  
  if (threshold) temp <- threshold(temp, tol)
  
  the_matrix <- drop0(temp)
  list(verified = nnzero(the_matrix) == 0, the_matrix=the_matrix) 
}

A_and_B_graph_seperated_by_S <- function(A, S, B, pcg, plot.subg=F) {
  # Test whether A and B are seperated in the graph (in PCG)
  all(B %in% all_seperated_from_A_by_S(A, S, pcg, plot.subg=plot.subg)$sep_set)
}

verifyL2MarkovPerfectness <- function(A, S, B, pcg, Sig, threshold=T, tol=DEFAULT.TOL, VERB=0, plot.subg=F) {
  g_verified <- A_and_B_graph_seperated_by_S(A, S, B, pcg, plot.subg = plot.subg)
  p_verified <- A_and_B_proj_seperated_by_S(A, S, B, Sig, threshold=threshold, tol=tol)$verified 
  pass_fail <- c("failed","passed")
  if (VERB) {
    print(sprintf("Graph %s, projection %s", pass_fail[g_verified+1], pass_fail[p_verified+1]))
  }
  
  (g_verified && p_verified) || (!g_verified && !p_verified)
    
}

verifyComponentDecomp <- function(j, S, pcg, Sig, all_nodes, restrict=F, plot_pcg=T, print.width=30, should_sort=F) {
  
  pcg_j_removed <- delete_vertices(pcg, j)
  #plot(pcg_j_removed, vertex.label=V(pcg_j_removed)$label)
  
  comps <-components(pcg_j_removed)
  G_comp <- list()
  for (k in 1:comps$no) {
    G_comp[[k]] <- V(pcg_j_removed)$label[comps$membership == k] 
  }
  
  S_comp <- lapply(G_comp, function(comp) intersect(comp,S) )
  
  cat(sprintf('%10s','Original:')) 
  orig_lat <- computeLattice(j, S, Sig, all_nodes, should_sort=should_sort)
  print( orig_lat, width=print.width, VERB=2) 
  cat('\n')
  Lat_comp <- list()
  for (k in 1:length(G_comp)){
    cat(sprintf('%10s',sprintf('S_%d:', k)))
    if (restrict) 
      Lat_comp[[k]] <- computeLattice(j, S_comp[[k]], Sig, all_nodes, 
                                      restricted_set=union(G_comp[[k]],j), should_sort=should_sort)  
    else
      Lat_comp[[k]] <- computeLattice(j, S_comp[[k]], Sig, all_nodes, should_sort=should_sort)  
    
    print(Lat_comp[[k]] , width=print.width, VERB=2) 
    cat('\n')
  }
  list_of_m <- lapply(Lat_comp, function(lat) lat$m)
  list_of_M <- lapply(Lat_comp, function(lat) lat$M)
  
  m_verified <- setequal( orig_lat$m, Reduce(union, list_of_m) )
  M_verified <- setequal( orig_lat$M, Reduce(union, list_of_M) )
  
  print.width.format <- sprintf('%%%ss',print.width)
  cat(sprintf('%10s','verified:'))
  cat(sprintf(print.width.format, if (m_verified) 'yes' else 'no'))
  cat(' -- ')
  cat(sprintf(print.width.format, if (M_verified) 'yes' else 'no'))
  cat('\n')
  #switch( isequal(Reduce(union,))   )
  
  if (plot_pcg) {
    plot(
      pcg_with_A_graphically_emphasized(
        pcg_with_S_graphically_removed(pcg,j)
      , S, "#CCCC80")
    )
    #labels <- V(pcg_j_removed)$label
    #colors <- rep("#FFFFAA",vcount(pcg_j_removed))
    #colors[labels %in% S] <- "#FFAAFF"
    #plot(pcg_j_removed,  vertex.label=labels, vertex.color=colors)  
  }
  
  list(verified = (m_verified && M_verified), 
       orig_lat=orig_lat, 
       Lat_comp=Lat_comp, 
       G_comp=G_comp, 
       S_comp=S_comp, 
       pcg_j_removed=pcg_j_removed)
  
}



#all_seperated_from_A_by_S <- function(A, S, pcg, all_nodes, plot_subg=F) {
all_seperated_from_A_by_S <- function(A, S, pcg, plot.subg=F) {
  # relies on pcg having node labels, i.e. V(pcg)$labes is defined
  
  #B <- setdiff(all_nodes, union(A,S))
  B <- setdiff(V(pcg)$label, union(A,S))
  
  pcg_subg <- delete_vertices(pcg,S)
  Dmat <-  distances(pcg_subg)
  labels <- V(pcg_subg)$label
  rownames(Dmat) <- labels
  colnames(Dmat) <- labels
  sub_Dmat <-  Dmat[labels %in% A, labels %in% B, drop=FALSE] 
  
  C <- B[colMins(sub_Dmat) == Inf]
  class(C) <- "set"
  if (plot.subg) {
    #plot(pcg_subg, vertex.label=V(pcg_subg)$labels)
    # temp_pcg <- pcg
    # labels <- V(temp_pcg)$label
    # E(temp_pcg)$lty <- 1
    # E(temp_pcg)$width <- 2
    # E(temp_pcg)[from(S)]$lty <- 2
    # E(temp_pcg)[from(S)]$width <- 0.5
    # colors <- rep("#FFFFAA",vcount(temp_pcg))
    # colors[labels %in% S] <- "#FFAAAA"
    # V(temp_pcg)$color <- colors
    # 
    # plot(temp_pcg)
    # #plot(temp_pcg, vertex.label=V(temp_pcg)$label, vertex.color=colors)
    plot(
      pcg_with_A_graphically_emphasized( 
        pcg_with_A_graphically_emphasized( 
          pcg_with_S_graphically_removed(pcg, S), 
          A, "#CCFFCC"),
        C, "#33FF33")
    )
    
  }
  
  list(sep_set=C, dist=sub_Dmat)
}

pcg_with_S_graphically_removed <- function(pcg,S) {
  labels <- V(pcg)$label
  E(pcg)$lty <- 1
  E(pcg)[from(S)]$lty <- 2
  E(pcg)$width <- 2
  E(pcg)[from(S)]$width <- 0.5
  
  V(pcg)$color[labels %in% S] <- "#FFAAAA"
  
  pcg
}

pcg_with_A_graphically_emphasized <- function(pcg, A, color) {
  labels <- V(pcg)$label
  V(pcg)$color[labels %in% A] <-  color #"#AAFFAA"
  pcg
}