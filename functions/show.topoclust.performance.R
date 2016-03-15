# This function receives a group of loci or gene trees. If necessary, it first estimates the trees, then it estimates the topological distance matrix. MDS is then used for the given number of dimesion values and the given number of ks for each. The output is a data matrix and a plot of how many dimensions were supported by each of the number of dimensions.

require(cluster)

topoclust.performance <- function(data, max.dim = 10, max.k = 10, makeplot = T, output.clusters = F){

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
	
	bestNclust <- vector()
	for(i in 2:max.dim){
		mdsres <- cmdscale(as.dist(topdistmat), eig = T, k = i)
		gapstats <- clusGap(mdsres$points, pam, max.k)
		bestNclust <- c(bestNclust, which(gapstats$Tab[, "gap"] == max(gapstats$Tab[, "gap"]))) 
		print(i)
	}
	
	resmat <- cbind(2:(length(bestNclust)+1), bestNclust)
	
	if(makeplot){
		plot(resmat, pch = 19, xlab = "Number of dimensions in MDS", ylab = "Number of clusters selected")
	}
	
	return(resmat)

} 