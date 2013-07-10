# TODO: Add comment
# 
# Author: mohamedahmed
###############################################################################


pathRanker <- function(gene.graph,K=10,minPathSize = 1,sameGenePenalty = 1/10,normalize = TRUE) {
	# Add S, T vetrices for the shortest path algorithm.
	gene.graph <- gene.graph + vertices("s", "t")
	gene.graph <- add.edges(gene.graph,
					c(rbind("s",V(gene.graph)[degree(gene.graph,mode="in")==1]$name), 
					rbind(V(gene.graph)[degree(gene.graph,mode="out")==1]$name,"t")), 
					attr=list(samegene=FALSE, 
							edge.weights=list(rep(-2, length(unlist(gene.graph$y.labels)) )), 
							attr=""))
	
	# Get edge.weights and apply ecdf on each column 
	edge.weights = do.call("rbind", E(gene.graph)$edge.weights)
	edge.probs <- apply(edge.weights, 2, function(x) ecdf(x[1:ecount(gene.graph)])(x)) 
	if (sum(E(gene.graph)$samegene) > 0) 
		edge.probs[E(gene.graph)$samegene,] <- sameGenePenalty
	
	# Normalize across Y labels
	if (ncol(edge.probs) > 1 & normalize == TRUE) {
		cat("Normalizing edges by response label.\n")
		edge.probs <- edge.probs / rowSums(edge.probs)
	}
	edge.probs <- -log(edge.probs)    
	
	# Prepare edgelist as required by PathRanker method.
	label = lapply(E(gene.graph)$attr, paste, collapse=' + ')	# For edges annotated with >1 compound, concatenate.
	edgelist = get.edgelist(gene.graph, names=FALSE)
	edgelist = data.frame(from=edgelist[,1], to=edgelist[,2], label=unlist(label), stringsAsFactors=FALSE)
	
	zret <- list()
	for (i in 1:ncol(edge.probs)) {
		if (ncol(edge.probs) > 1) cat("Extracting the",K,"most probable paths for",gene.graph$y.labels[i],"\n")
		else cat("Extracting the",K,"most probable paths.\n")
		
		
		# min path size is increased by 2 to include start and end compounds!!!
		ps <- .Call("pathranker",
				V(gene.graph)$name, 	#nodes
				edgelist,				#Edges
				edge.probs[,i],			#Edge weights
				K=K,minpathsize = minPathSize+2)
		
		
		idx <- which(sapply(ps,is.null))		
		if (length(idx)>0) ps <- ps[-idx]
		
		if (ncol(edge.probs) > 1) zret[[i]] <- ps
		else zret <- ps
	}
	
	if (ncol(edge.probs) > 1) {
		names(zret) <- gene.graph$y.labels
		colnames(edge.probs) <- paste("prob",gene.graph$y.labels,sep = ":")
	} else names(edge.probs) <- "prob"
	
	return(list(paths = zret,edge.weights = edge.probs))
}

