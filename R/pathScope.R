samplePaths <- function(rnet,edge.weights,maxpathlength,numberofsamples,numberofwarmup) {
    if (maxpathlength > nrow(rnet$edges)) {
        cat("\nMaximum Path Length > Number of Network Edges ... setting them to be equal")
        maxpathlength <-  nrow(rnet$edges)
    }


    gedges <- data.frame(from = as.integer(rnet$edges$from), 
                         to = as.integer(rnet$edges$to),
                         label =  paste(rnet$nodes[rnet$edges$from],rnet$nodes[rnet$edges$to],sep = ":"),
                         pathway = rep(rnet$sbml.info$network,nrow(rnet$edges)),
                         stringsAsFactors = FALSE)

    ps <- .Call("samplepaths",
                rnet$nodes,
                gedges,
                edge.weights,
                MAXPATHLENGTH = as.integer(maxpathlength),
                SAMPLEPATHS = as.integer(numberofsamples),
                WARMUPSTEPS = as.integer(numberofwarmup))
    return(matrix(ps,maxpathlength+1,numberofsamples,byrow = TRUE))
}

scope <- function(geneNetwork,edge.weights,sampledpaths,alpha = 0.01,echo = 0) {    
    z <- .Call("scope",
            NODES = geneNetwork$nodes,
            EDGES = geneNetwork$edges,
            WEIGHTS = -log(edge.weights),
            SAMPLEDPATHS = t(sampledpaths),
            ALPHA = alpha,
            ECHO = as.integer(echo))
    zidx <- which(sapply(z$paths,is.null) == FALSE)
    z$paths <- z$paths[zidx]
    return(z)
}
