# This gets the base composition for every taxon for a group of loci. The result is a list of matrices of base compositions.

get.taxa.base.comps <- function(loci){

	cont <- list()
	for(i in 1:length(loci)){
	      cont[[i]] <- matrix(NA, nrow = nrow(birdat[[1]]), ncol = 4)
	      for(j in 1:nrow(loci[[i]])){
	      	    cont[[i]][j,] <- table(as.character(loci[[i]])[j,])[c("a", "c", "g", "t")]
		    cont[[i]][j,] <- cont[[i]][j,] / sum(cont[[i]][j,])
	      }
	      rownames(cont[[i]]) <- rownames(loci[[i]])
	      print(paste("processed locus", i))
	}

	return(cont)
	
}