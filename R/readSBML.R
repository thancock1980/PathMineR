readSBML <- function(filename) {
    if (!file.exists(filename)) stop("Cannot find file!")

    zsbml <- .Call("readsbmlfile", FILENAME = filename)

    edge.list <- do.call("rbind",lapply(lapply(zsbml$reaction.info,"[[","edges"),"data.frame",stringsAsFactors = FALSE))
    dimnames(edge.list)[[1]] <- 1:nrow(edge.list)
    
	reaction.id <- names(zsbml$reaction.info)
	galines <- read.delim("ga.txt", header=FALSE, stringsAsFactors=FALSE, skip=1)
	galines <- galines[,c(1,2,3,5,6,16)]
	# Reformatting the list, removing unnecessary info.
	galines[,5] <- gsub("REACTOME:","", galines[,5])
	galines[,6] <- gsub("(.+):(\\w+)(\\)+)","\\2",galines[,6], perl=TRUE)
	names(galines) <- c("Type", "ID", "Name", "GO", "Reaction", "Complex")
	
	
	reaction.info <- lapply(reaction.id,function(x) unique(galines[galines[,5]==x & galines[,1]=='UniProtKB' ,2]))
	names(reaction.info) <- reaction.id
	unique.genes <- unique(unlist(reaction.info))
	names(unique.genes) <- unique.genes
	gene.list = lapply(unique.genes, function (x) data.frame(galines[galines[,2]==x & galines[,5] %in% reaction.id,c(-1,-2)], row.names=NULL))
	
	
	#reaction.info <- lapply(reaction.id,function(x) droplevels(galines[galines[,5]==x & galines[,1]=='UniProtKB' ,2]))
	
	#lapply(reaction.names, function(x) {z$sbml.info$reaction.info[[x]]$genes <<- reaction.info[[x]]} )
	#zsbml$sbml.info$reaction.info[]
	
	zsbml$species.info <- data.frame(zsbml$species.info,stringsAsFactors = FALSE)
	#print("Im here")
	
    kobj <- list(edge.list = edge.list,
				 reaction.info = reaction.info,
				 gene.list = gene.list,
				 species.info = zsbml$species.info,
				 sbml.info = zsbml)
   
   return(kobj)
}

makeIgraph <- function(reaction.substrate.network){
	nodes <- reaction.substrate.network$nodes
	edges <- reaction.substrate.network$edges
	z <- reaction.substrate.network$z
	
	# contructing a bipartite graph
	graph = graph.bipartite( grepl("REACT_",nodes), c(t(as.matrix(edges)-1)), directed=TRUE)
	V(graph)$shape<- ifelse(V(graph)$type==TRUE, "square", "circle")
	V(graph)$color <- ifelse(V(graph)$type==TRUE,"red", "skyblue")
	V(graph)$size <- 5
	E(graph)$arrow.size = 0.2
	
	# Exporting Species & reactions attributes
	reactions = sapply(nodes[grepl("REACT_", nodes)], function(x) z$sbml.info$reaction.info[x], USE.NAMES=FALSE)
	
	names(nodes) <- nodes
	species = lapply(nodes[grepl("spec", nodes)], function(x)as.list(z$sbml.info$species.info[z$sbml.info$species.info[,1]==x,]))
	
	attr = c(reactions,species)
	lapply(names(attr), function(x) V(graph)[match(x, net$nodes)-1]$attr <<- attr[x])
	graph$sbml = z
	return(graph)
}


makeReactionSubstrateNetwork <- function(kobj){
	
	reactants <- lapply(kobj$sbml.info$reaction.info, "[[", "reactants")
	products <- lapply(kobj$sbml.info$reaction.info, "[[", "products")
	
	edge.list <- do.call("rbind",lapply(names(reactants), 
					function(x) expand.grid(reactants[[x]], x, stringsAsFactors=FALSE)))
	edge.list <- rbind(
					do.call("rbind",lapply(names(products), 
					function(x) expand.grid(x, products[[x]], stringsAsFactors=FALSE)))
				 , edge.list)
	

	names(edge.list) <- c("from","to")
	# edges with the same start and end
	rmidx <- which(edge.list$from == edge.list$to)
	if (length(rmidx) > 0) edge.list <- edge.list[-rmidx,]
	
	# Remove reactions with only products or reactions
	# These reactions are annotated to be the expression of vertain proteins and/or modifiers
	# thus irrelevant to the reaction-substrate network
	rmidx <- which(edge.list$from %in% setdiff(names(products), edge.list$to))
	edge.list <- edge.list[-rmidx,]
	
	#make the network
	gnds <- unique(unlist(edge.list, use.names=FALSE))
	edge.list$from  <- match(edge.list$from,gnds)
	edge.list$to  <- match(edge.list$to,gnds)
	
	return(list(nodes = gnds,edges = edge.list, z = kobj))	
}


makeGeneNetwork2 <- function(reaction.network){
	reaction.info <- reaction.network$reaction.info
	nodes <- reaction.network$nodes
	edges <- reaction.network$edges
	edge.list <- data.frame()
	
	for(i in 1:length(nodes)){
		from <- edges[edges$from == i, 2]
		list = lapply(from, function(x) expand.grid(reaction.info[[nodes[i]]], reaction.info[[nodes[[x]]]], stringsAsFactors=FALSE))
		
		edge.list <- rbind(do.call("rbind", lapply(list , "data.frame", stringsAsFactors=FALSE)), edge.list)
	}
	
	names(edge.list) <- c("from","to")
	# edges with the same start and end
    rmidx <- which(edge.list$from == edge.list$to)
    if (length(rmidx) > 0) edge.list <- edge.list[-rmidx,]
	
	#make the network
    gnds <- unique(unlist(reaction.info))
    edge.list$from  <- match(edge.list$from,gnds)
    edge.list$to  <- match(edge.list$to,gnds)
 	
 	return(list(nodes = gnds,edges = edge.list,reaction.info = reaction.info,sbml.info = reaction.network$sbml.info))   
}

makeReactionNetwork2 <- function(kobj,removeSmallCompounds = NULL) {
    redges <- kobj$edge.list
    edges <- NULL
    
    reaction.info <- kobj$reaction.info
    rlens <- sapply(reaction.info,length)
    nullreact <- names(reaction.info)[which(rlens == 0)]

    #remove reactions with no genes
    rmidx <- which(redges$label %in% nullreact)
    if (length(rmidx) > 0) redges <- redges[-rmidx,]
    all.reactions <- unique(redges$label)

    reaction.info <- reaction.info[all.reactions]

    rsbml <- kobj$sbml.info
    rsbml$reaction.info <- rsbml$reaction.info[all.reactions]

    # remove the irrelevant edges (if any)
    if (!is.null(removeSmallCompounds)) {
      rmidx <- c()
      for (i in 1:length(removeSmallCompounds)) {      
        rmidx <- sort(unique(c(rmidx,grep(removeSmallCompounds[i],redges$from),grep(removeSmallCompounds[i],redges$to))))
      }      
      if (length(rmidx) > 0) {
        redges <- redges[-rmidx,]
      }
    }

    for (reaction in all.reactions) {
        # all compounds (to.comps) produced by reaction
        to.comps <- unique(redges[redges$label == reaction,]$to)        

        # all reactions starting from compounds (to.comps)
        r.to <- NULL
        if (!is.null(to.comps)) r.to <- unique(redges[redges$from %in% to.comps,]$label)
    
        # if two reactions have different names but connect the same compounds in the same way with the same genes
        # then remove them
        if (length(r.to) > 0) {
            indx <- 1; rmidx <- c()
            for (r in r.to) {
                etest <- identical(redges[redges$label == reaction,][1:2],
                                   unique(rbind(redges[redges$label == reaction,][1:2],redges[redges$label == r,][1:2]))) 
                gtest <- identical(reaction.info[[reaction]],unique(c(reaction.info[[reaction]],reaction.info[[r]])))
                if (etest & gtest) rmidx <- c(rmidx,indx)
                indx <- indx + 1
            }
            if (length(rmidx) > 0) r.to <- r.to[-rmidx]
        }

        # if an edge exists then...
        rowframe <- NULL
        if (length(r.to) > 0) rowframe <- data.frame(from = reaction,to = r.to,stringsAsFactors = FALSE)
       
        if (!is.null(rowframe)) {    
            if (is.null(edges)) edges <- rowframe
            else edges <- rbind(edges,rowframe)
        }
    }
    # incase edges are duplicted
    edges <- unique(edges)

    # edges with the same start and end
    rmidx <- which(edges$from == edges$to)
    if (length(rmidx) > 0) edges <- edges[-rmidx,]

    edges <- data.frame(edges,stringsAsFactors = FALSE)
    names(edges) <- c("from","to")

    #make the network
    gnds <- names(reaction.info) 
    edges$from  <- match(edges$from,gnds)
    edges$to  <- match(edges$to,gnds)

    reaction.info <- lapply(reaction.info,toupper)
   
    return(list(nodes = gnds,edges = edges,reaction.info = reaction.info,sbml.info = rsbml))
}

computeMaxReactionCorrelations <- function(reaction.network,microarray,y=NULL) {
    zedgelist <- reaction.network$edges
    zedgelist$from <- reaction.network$nodes[zedgelist$from]
    zedgelist$to <- reaction.network$nodes[zedgelist$to]

    # check for reverse edges
    zforward <- paste(zedgelist$from,zedgelist$to,sep = ":")
    zreverse <- paste(zedgelist$to,zedgelist$from,sep =":")
    
    # make sure the first repeat is computed first, and the second just looked up
    zrepeat <- as.numeric(sapply(zforward,match,zreverse)) 
    gm <- which(1:length(zrepeat) < zrepeat)
    if (length(gm) > 0) zrepeat[gm] <- -1
    zrepeat[is.na(zrepeat)] <- -1
    zedgelist$reverseindex <- as.integer(zrepeat-1)

    # convert edge list to reaction indices
    zreactionlist <- reaction.network$reaction.info

    zedgelist$from <- as.integer(match(zedgelist$from,names(zreactionlist))-1)
    zedgelist$to <- as.integer(match(zedgelist$to,names(zreactionlist))-1)

    # convert reaction gene name list to indices within the microarray

    zreactiongeneindicies <- lapply(lapply(lapply(zreactionlist,match,names(microarray)),"-",1),as.integer)

    cors <- NULL
    if (!is.null(y)) {
        cors <- list()
        for (i in 1:nlevels(y)) {
            cat("\nComputering Correlations for ",levels(y)[i])
            cors[[i]] <- .Call("compute_max_reaction_correlation",zedgelist,zreactiongeneindicies,as.matrix(microarray[y == levels(y)[i],]))
        }
        cors <- data.frame(do.call("cbind",cors))
        names(cors) <- paste("cor",levels(y),sep = ":")
    }
    else {
        cat("\nComputering Correlations")
        cors <- .Call("compute_max_reaction_correlation",zedgelist,zreactiongeneindicies,as.matrix(microarray))
        cors <- data.frame(cors)
        names(cors) <- "cor"
    }
    cat("\n")

    return(cors)
}

defineReactionNeighborhood <- function(rnet,edge.weights,seed,range,ann = "GRFT",mode = "out") {
    sreactions <- which(sapply(lapply(lapply(rnet$sbml.info$reaction.info,"[[","reactants"),"%in%",seed),sum)>0)
    if (length(sreactions) == 0) {
      cat("\nSeed Compound",seed,"is never a substrate")
       return(NULL)
    } 

    # define start node sets   
    sreactions <- names(rnet$sbml.info$reaction.info)[sreactions]
    sridx <- match(sreactions,rnet$nodes)
    sedges <- data.frame(from = rep(1,length(sreactions)),to = sridx+1)
    gedges <- rbind(sedges,rnet$edges+1)    
    nds <- c("s",rnet$nodes)
    sweights <- data.frame(matrix(1,nrow = nrow(sedges),ncol = ncol(edge.weights)))
    names(sweights) <- names(edge.weights)
    weights <- rbind(sweights,edge.weights)

    # get the neighborhood from seed node
    ge <- graph.edgelist(el = as.matrix(gedges-1),directed = TRUE)

    nidx <- neighborhood(ge,order = range,nodes = V(ge)[0],mode = mode)[[1]]
    if (sum(nidx==0) > 0) nidx <- nidx[-which(nidx==0)]
    if (length(nidx) == 0) {
      cat("\nNeighborhood for",seed,"for mode",mode,"and range",range,"not valid")
      return(NULL)
    }
    scp <- nds[nidx]
       
    # define end node sets
    tedges <- data.frame(from = nidx+1,to = rep(length(nds)+1,length(nidx)))
    gedges <- rbind(gedges,tedges) 
    nds <- c(nds,"t")
    tweights <- data.frame(matrix(1,nrow = nrow(tedges),ncol = ncol(edge.weights)))
    names(tweights) <- names(edge.weights)
    weights <- rbind(weights,tweights)
    
    glabel <- paste(nds[gedges$from],nds[gedges$to],sep = ":")
    gpath <- rep(rnet$sbml.info$network,nrow(gedges))
    gedges$from <- as.integer(gedges$from)
    gedges$to <- as.integer(gedges$to)
    gedges$label <- glabel
    gedges$pathway <- gpath
    
    nobj <- list(nodes = nds,edges = gedges)  

    return(list(geneNetwork = nobj,edge.weights = weights,seed = seed,neighborhood = scp))
}



