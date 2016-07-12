# This function takes the 100-tip tree from the mcmctree output (with confidence intervals) and calculates the accuracy and uncertainty of two nodes with known age: the root (age = 100), and N2, parent to t1 to 10 (age = 10).

get.mcmctree.acc.unc <- function(phy){
	brtimes <- branching.times(phy)
	N2 <- getMRCA(phy, paste0("t", 1:10)) - 100
	N2age <- brtimes[N2]
	N2acc <- (N2age - 0.1) / 0.1

	N2CI <- diff(as.numeric(strsplit(phy$node.label[N2], "-")[[1]]))
	N2unc <- N2CI/0.1

	rootage <- brtimes[1]
	rootacc <- (rootage - 1) / 1

	rootCI <- diff(as.numeric(strsplit(phy$node.label[1], "-")[[1]]))
	rootunc <- rootCI/1

	output <- matrix(c(N2acc, rootacc, N2unc, rootunc), 2, 2, dimnames = list(c("N2", "root"), c("Accuracy", "Uncertainty")))
	
	return(output)
}