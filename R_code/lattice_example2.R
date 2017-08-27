### Source the base code for computing SEM
#source("gaussian_projections_base.R") # make sure this file is in your current directory, or add /path/to/ as needed

source('load_lattice_data3.R')

# check whether PCG mathces Gamma
all ( (abs(Gamma) > 0)*1 == (get.adjacency(pcg) + diag(1,d)) )

range(Sig)
range(Gamma)

#image(Bet, sub="", xlab="", ylab="")
image(Gamma, sub="", xlab="", ylab="")
image(Sig, sub="", xlab="", ylab="")

# Plot DAG
#plot(dag)

# Plot PCG
plot(pcg)

print( getParentCoefs2(2, all_nodes[-2], Sig, threshold=T)$coefs,  col.names=T)

j <- 4
S <- all_nodes[-j]
lat <- computeLattice(j, c(3,4,10), Sig, all_nodes, should_sort=F)
lat
isInLattice(c(6,3,8,4,10),lat)
print( computeLattice(j, c(1,8,9,5,7), Sig, all_nodes), width=30)


lattices <- computeAllSubordinateLattices(j, Sig, all_nodes, print.width = 30, VERB=2)
#is_set_in_lattices(c(1,3,7,8,9),lattices)
S <- c(2,3,4)
is_set_in_lattices(S,lattices)
lattices <- computeAllSubordinateLattices(j, Sig, all_nodes, root_S = S, print.width = 30, VERB=2)


S <- sample(1:d, sample(d-1,1))
j <- sample(setdiff(all_nodes, S),1)

S <- c(7,3,8,10)
j  <- 5
out <- verifyComponentDecomp(j, S, pcg, Sig, all_nodes, should_sort=T)


S <- c()
A <- 2
B <- c(1,6,9,12,16)
A_and_B_graph_seperated_by_S(A, S, B, pcg, plot.subg=T)
A_and_B_proj_seperated_by_S(A, S, B, Sig, threshold=1) 
verifyL2MarkovPerfectness(A, S, B, pcg, Sig, VERB=1)


A <- c(10)
j <- 2
R <- c(1,3,5,7,9)
verifyL2MarkovPerfectness(j, A, R, pcg, Sig, VERB=1, plot.subg=T)


# Problematic ones
j <- 1
S <- sample(2:d, sample(d-1,1))
out <- verifyComponentDecomp(j, S, pcg, Sig, all_nodes, should_sort=T)

