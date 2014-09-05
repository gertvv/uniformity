# Given a matrix of coordinate vectors (rows), calculate the Minimum Spanning
# Tree.
# Implemented for memory efficiency: does NOT store distance matrix.
# So memory use is O(n), not O(n^2). Time complexity O(n^2).
euclidMST <- function(mat) {
	.C("uniformity_mst", as.double(mat), as.integer(dim(mat)[1]), as.integer(dim(mat)[2]), result=integer(dim(mat)[1] - 1))$result
}

## Test that "euclidMST" is correctly implemented:
# library(vegan)
# all(euclidMST(bench) == spantree(dist(bench))$kid)

# Given a matrix of coordinate vectors (rows), calculate the nearest
# neighbours.
# Implemented for memory efficiency: does NOT store distance matrix.
# So memory use is O(n), not O(n^2). Time complexity O(n^2).
euclidNN <- function(mat) {
	res <- .C("uniformity_nn", as.double(mat), as.integer(dim(mat)[1]), as.integer(dim(mat)[2]), nearest=integer(dim(mat)[1]), sqdist=double(dim(mat)[1]))
	list(nearest=res$nearest, distance=sqrt(res$sqdist))
}

## Test that "euclidNN" is correctly implemented:
# shortest1 <- apply(as.matrix(dist(bench)) + diag(NA, dim(bench)[1]), 1, function(x) { min(x, na.rm=T) })
# shortest2 <- euclidNN(bench)$distance
# all(shortest1 == shortest2)

# Friedman-Rafsky test that two samples are from the same population (null).
# [Friedman and Rafsky, 1979 and Smith and Jain, 1984]
testUniformMST <- function(bench, test) {
	# Generate connectivity list: if st[i] = j, there is an edge (i + 1) -- (j)
	st <- euclidMST(rbind(bench, test))
	# Transform the list to obtain true edge-pairs
	st <- t(sapply(1:length(st), function(i) { c(i + 1, st[i]) }))

	N <- dim(bench)[1]
	M <- dim(test)[1]
	L <- M + N
	# Expected value of T (cross-connectivity)
	E <- 2*M*N/L
	# Calculate cross-connectivity between the sets
	T <- sum(apply(st, 1, function(e) { if ((e[1] <= N && e[2] > N) || (e[1] > N && e[2] <= N)) { 1 } else { 0 } }))

	# C is the number of edge pairs that share a vertex. Calculating it directly (like this) is slow:
	## C <- sum(sapply(1:(dim(st)[1] - 1), function(i) { e <- st[i,]; rem <- tail(st, -i); sum(rem[,1] == e[1] | rem[,1] == e[2] | rem[,2] == e[1] | rem[,2] == e[2]) }))
	# But it is related to the degree of the vertices, determined as:
	d <- table(c(st[,1], st[,2]))
	C <- sum(d * (d - 1) / 2)

	# Variance of T
	V <- (2 * M * N / (L * (L - 1))) * ((2 * M * N - L) / L + (C - L + 2) / ((L - 2) * (L - 3)) * (L * (L - 1) - 4 * M * N + 2))

	# values that can be then Z-tested(?)
	(T - E)/sqrt(V)
}

# Test uniformity using the Coefficient of Variation of the nearest-neighbour
# distances
testUniformCOV <- function(data) {
	# Calculate minimum inter-vertex distances
	gamma <- euclidNN(data)$distance

	# Calculate COV measure
	sd(gamma)/mean(gamma)
}

# Test uniformity using the Mesh Ratio
testUniformMR <- function(data) {
	# Calculate minimum inter-vertex distances
	gamma <- euclidNN(data)$distance

	# Calculate MR measure
	max(gamma) / min(gamma)
}

# Test uniformity by normalized error in mean (centroid)
testUniformComponentErrors <- function(data, centroid) {
	abs(apply(data, 2, mean) - centroid) / apply(data, 2, sd)
}

testUniformRatioError <- function(data, centroid) {
	n <- length(centroid)
	rr <- centroid[1] / centroid[n] # The real ratio
	est <- apply(data, 2, mean)
	er <- est[1] / est[n] # The estimated ratio
	abs(rr - er) / rr # Relative error
}
