# Process a pair of simulated chronograms. Take pair of chronogs in a list. Simulate data from them. Re-estimate phylograms. Output gene tree list.

require(phangorn)
require(NELSI)

process.gd.sim <- function(trees, nclocks = 3, meanclocks = c(0.01, 0.015, 0.02), sdclocks = c(0.01, 0.01, 0.01), seqlength = 100, reps = 100){

	       # Simulate phylograms
	       if(length(trees) * nclocks * reps != 200) reps <- 200 / (length(trees) * nclocks)
	       allsims <- list()
	       for(i in 1:nclocks){
	       	     for(j in 1:length(trees)){
		     	   allsims[[length(allsims) + 1]] <- list()
			   names(allsims)[length(allsims)] <- paste0("tr", j, "pm", i)
			   trees[[j]]$edge.length <- rlnorm(length(trees[[j]]$edge.length), meanlog = log(meanclocks[i]), sdlog = sdclocks[i])
			   for(k in 1:reps){
			   	 allsims[[length(allsims)]][[k]] <- trees[[j]]
				 scaleval <- rlnorm(1, log(1), 0.25)
				 allsims[[length(allsims)]][[k]]$edge.length <- allsims[[length(allsims)]][[k]]$edge.length * scaleval
			   }
		     }
	       }
	       
	       # Simulate alignment data from phylograms, and estimate trees.
	       
	       simals <- list()
	       estrees <- list()
	       for(i in 1:length(allsims)){
	       	     print(names(allsims[i]))
		     simals[[i]] <- lapply(allsims[[i]], simSeq, l = seqlength)
		     estrees[[i]] <- lapply(simals[[i]], function(x) optim.pml(pml(NJ(dist.dna(as.DNAbin(x))), x), optNni = T)$tree)
		     
		     print(paste("###############################  COMPLETED SCHEME", i,",", names(allsims)[i],"!!!!!! ####################################"))
	       }

	       names(simals) <- names(allsims)
	       names(estrees) <- names(allsims)
	       
	       reslist <- list(sim.phylogs = allsims, sim.genes = simals, est.phylogs = estrees)

	       return(reslist)
}