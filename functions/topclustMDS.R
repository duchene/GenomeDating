# This function receives a group of loci or gene trees. If necessary, it first estimates the trees, then it estimates the topological distance matrix. MDS is then used for the given number of dimesion values and the given number of ks for each. The output is a data matrix and a plot of how many dimensions were supported by each of the number of dimensions.

require(cluster)

topoclustMDS <- function(data, mdsdim = 1, makeplot = T){
	
	max.k = length(data) - 1
	
	if(class(data[[1]]) == "phyDat"){
		als <- data
		data <- lapply(als, function(x) optim.pml(pml(NJ(dist.dna(as.DNAbin(x))), x), optNni = T)$tree)
	}
	topdistmat <- matrix(NA, ncol = length(data), nrow = length(data))
	for(i in 1:length(data)){
	      for(j in i:length(data)){
	      	    topdistmat[j, i] <- dist.topo(data[[i]], data[[j]])
	      }
	}
	
	print("Created topological distance matrix")
	#print(topdistmat)
	
	mdsres <- cmdscale(as.dist(topdistmat), eig = T, k = mdsdim)
	#print(head(mdsres$points))
	print("MDS completed")
	if(makeplot) plot(mdsres$points)
	gapstats <- clusGap(mdsres$points, pam, max.k)
        bestNclust <- which(gapstats$Tab[, "gap"] == max(gapstats$Tab[, "gap"]))
	print("Gap statistics table:")
	print(gapstats$Tab)
	if(length(bestNclust) > 1) bestNclust <- bestNclust[1]
	clusterdata <- pam(mdsres$points, k = bestNclust)
	
	res <- list(topodists = topdistmat, mds = mdsres, clustering.data = clusterdata, k = bestNclust)
	
	return(res)
	
} 