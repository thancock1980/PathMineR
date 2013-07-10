# TODO: Add comment
# 
# Author: ahmedmohamed
###############################################################################
KGML2igraph <- function(filename){
	zkgml = .Call("readkgmlfile", FILENAME = filename)
	
}

SBML2igraph <- function(filename){
	if (!file.exists(filename)) stop("Cannot find file!")
	
	zsbml <- .Call("readsbmlfile", FILENAME = filename)
	
	print("SBML file processed successfully :|")
	print("Processing Gene Association file: Slow,will be removed in the final version")
	#################Gene Association File ###################
	ga.file = paste(find.package("MetPathMineR"),"/extdata/ga.txt", sep="")
	galines <- read.delim(ga.file, header=FALSE, stringsAsFactors=FALSE, skip=1)
	galines = galines[galines[,13]=="taxon:9606",] #Restrict to Homo sapiens
	galines <- galines[,c(1,2,3,5,6,16)]
	# Reformatting the list, removing unnecessary info.
	galines[,5] <- gsub("REACTOME:","", galines[,5])
	galines[,6] <- gsub("(.+):(\\w+)(\\)+)","\\2",galines[,6], perl=TRUE)
	names(galines) <- c("Type", "ID", "Name", "GO", "Reaction", "Complex")
	##########################################################
	
	print("File Processed")
	## Manipulating genes in the reaction network#############
	species.info <- data.frame(zsbml$species.info,stringsAsFactors = FALSE)
	reaction.info <- zsbml$reaction.info
	
	reaction.id <- names(reaction.info)
	lapply(reaction.id,function(x) reaction.info[[x]]$genes <<-		#Combine genes from SBML and Reactome list files.
								      unique(c(species.info[species.info$id %in% unlist(reaction.info[[x]]$genes),8],
											  galines[galines[,5]==x & galines[,1]=='UniProtKB' ,2]
		  								)))
	genes = unique(unlist(lapply(reaction.info, "[[", "genes")))
	names(genes) <- genes
	
	gene.list = lapply(genes, function (x) data.frame(galines[galines[,2]==x & galines[,5] %in% reaction.id,c(-1,-2)], row.names=NULL))
	##########################################################
	
		
	print("Constructing Reaction Substrate Network")
	######### Make the reaction substrate network ############ 
	reactants <- lapply(reaction.info, "[[", "reactants")
	products <- lapply(reaction.info, "[[", "products")
	
	edges <- do.call("rbind",lapply(names(reactants), 
					function(x) expand.grid(reactants[[x]], x, stringsAsFactors=FALSE)))
	edges <- rbind(
			do.call("rbind",lapply(names(products), 
							function(x) expand.grid(x, products[[x]], stringsAsFactors=FALSE)))
			, edges)
		
	names(edges) <- c("from","to")
	# edges with the same start and end
	rmidx <- which(edges$from == edges$to)
	if (length(rmidx) > 0) edges <- edges[-rmidx,]
	
	#Remove small compounds  #Reactome only 
	small.compounds.list <- readLines("small_comps.list")
	small.compound.ids <- species.info[species.info[,2] %in% small.compounds.list, 1]
	edges <- edges[ !(edges$from%in%small.compound.ids | edges$to%in%small.compound.ids),]
	
	# Remove reactions with only products or reactions
	# These reactions are annotated to be the expression of vertain proteins and/or modifiers
	# thus irrelevant to the reaction-substrate network
	rmidx <- which(edges$from %in% setdiff(names(products), edges$to))
	if (length(rmidx) > 0) edges <- edges[-rmidx,]
	
#	#make the network
#	nodes <- unique(unlist(edges, use.names=FALSE))
#	edges$from  <- match(edges$from,nodes)
#	edges$to  <- match(edges$to,nodes)
#	#######################################################
	
	########## Making the iGraph object ###################
	# contructing a bipartite graph
	graph = graph.data.frame(edges)
	V(graph)$type <- grepl("REACT_",V(graph)$name)
	V(graph)$shape<- ifelse(V(graph)$type==TRUE, "square", "circle")
	V(graph)$color <- ifelse(V(graph)$type==TRUE,"red", "skyblue")
	V(graph)$size <- 5
	E(graph)$arrow.size <- 0.2
	#V(graph)$id <- nodes	
	#nodes <- V(graph)$name
		

	# Exporting Species & reactions attributes
	#reactions <- sapply(nodes[grepl("REACT_", nodes)], function(x) reaction.info[x], USE.NAMES=FALSE)
#	V(graph)[type]$attr <- reaction.info[V(graph)[type]]	
	#names(nodes) <- nodes	
	
	#species <- lapply(nodes[grepl("spec", nodes)], function(x)as.list(species.info[species.info[,1]==x,]))
	#V(graph)[!type]$attr <- species[V(graph)]
	print("Assigning node attributes: Slow, for speed, C function needs to return a named list instead of data.frame")
		
	reactions <- sapply(V(graph)[type]$name, function(x) reaction.info[x], USE.NAMES=FALSE)
	V(graph)[names(reactions)]$attr <- reactions
	V(graph)[type]$genes = lapply(V(graph)[type]$attr, "[[", "genes") 
	
	species <- lapply(V(graph)[!type]$name, function(x) as.list(species.info[species.info[,1]==x,]))
	names(species) = lapply(species, "[[", "id")
	V(graph)[names(species)]$attr <- species

	
	
	attr <- c(reactions,species)
	lapply(names(attr), function(x) V(graph)[x]$attr <<- attr[x])
	#######################################################
	
	graph$gene.list <- gene.list
	graph$attr <- attr
	graph$reaction.info <- reaction.info
	graph$species.info <- species.info
	#graph$reactions = reactions
		
	return(graph)
}

makeReactionNetwork <- function(graph, simplify=FALSE){
	edges <- get.edgelist(graph)

	
	input <- edges[edges[,2]%in% V(graph)[type]$name,]  #Reactant-reaction edges
	output <- edges[edges[,1]%in% V(graph)[type]$name,]   #Reaction-product edges
	
	match1 <- match(output[,2],input[,1])    		#Find reactions that produces reactants in input 
	match2 <- match(input[,1],output[,2])			#Find reactions that consumes products in output
	
	#Connect matches
	reaction.connections <- rbind(						
		cbind(output[,],input[match1,2]),
		cbind(output[match2,], input[,2])
	)
	
	
	
	# remove missing and duplicate values 
	reaction.connections <- reaction.connections[complete.cases(reaction.connections),]
	reaction.connections <- unique(reaction.connections)
	
	
	reaction.edges <- reaction.connections[,c(1,3)]
	reaction.edges <- unique(reaction.edges)
	rmidx <- which(reaction.edges[,1]==reaction.edges[,2])			#Remove self-edges
	if (length(rmidx) > 0) reaction.edges <- reaction.edges[-rmidx,]
	
	
	#### Making a reaction network graph from edges#####
	reaction.graph <- graph.edgelist(reaction.edges)
#	V(reaction.graph)$id <- V(graph)[type==TRUE]$id
	V(reaction.graph)$attr <- V(graph)[V(reaction.graph)$name]$attr
	V(reaction.graph)$genes <- V(graph)[V(reaction.graph)$name]$genes
	V(reaction.graph)$size <- 5
	E(reaction.graph)$arrow.size <- 0.2
	E(reaction.graph)$attr <- apply(reaction.edges, 1, 
			function(x) reaction.connections[reaction.connections[,1]==x[1] & reaction.connections[,3]==x[2], 2])
	
	reaction.graph$species.info <- do.call("rbind", lapply(V(graph)[type==FALSE]$attr, "data.frame", stringsAsFactors=FALSE))	
	reaction.graph$gene.list <- graph$gene.list	
	reaction.graph$reaction.coor = cbind(reaction=1:vcount(reaction.graph), 
							met.reaction=match(V(reaction.graph)$name, V(graph)$name) )
#	reaction.graph$reaction.connections <- reaction.connections
	if(simplify)
		reaction.graph <- simplifyReactionNetwork(reaction.graph, missing.genes=FALSE)
	return(reaction.graph)
}


makeGeneNetwork <- function(reaction.graph){
	edges <- get.edgelist(reaction.graph)
	genes <- V(reaction.graph)$genes
	names(genes) <- V(reaction.graph)$name
	genes = mapply(paste, genes, names(genes), sep="#") #Duplicate Gene edges
	
	# Convert reaction edges to gene edges.
	gene.connections <- do.call("rbind" ,
			lapply(1:dim(edges)[1], function(x) 
							expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)
	))
	colnames(gene.connections)	<- c("from", "to", "id")

	#make the graph.
	gene.graph = graph.data.frame(gene.connections[,-3])
	V(gene.graph)$attr <- reaction.graph$gene.list[sub("#(.+)$", "",V(gene.graph)$name)]
	V(gene.graph)$size <- 5
	E(gene.graph)$arrow.size <- 0.2
	
	# get the edge attributes
	E(gene.graph)$attr <- E(reaction.graph)$attr[gene.connections$id]
	E(gene.graph)$samegene <- sub("#(.+)$", "",gene.connections$from) == sub("#(.+)$", "",gene.connections$to)  
	
	# Map vertex indices to the reaction graph (for the layout function). 
	gene.coor <- cbind(gene=1:vcount(gene.graph), 
					reaction=match(sub("^(.+)#", "",V(gene.graph)$name), V(reaction.graph)$name))
	
	# get graph attributes.
	gene.graph$species.info <- reaction.graph$species.info
	reaction.info <- V(reaction.graph)$attr
	names(reaction.info) <- V(reaction.graph)$name
	gene.graph$reaction.info <- reaction.info
	
	return(gene.graph)
}

simplifyReactionNetwork <- function(reaction.graph, missing.genes=TRUE, reconnect.threshold=vcount(reaction.graph)){
	#List of all missing gene reactions.
	V(reaction.graph)$no.gene = lapply( V(reaction.graph)$genes, length) == FALSE
	species.info = reaction.graph$species.info
	translocation = sapply( V(reaction.graph)[no.gene]$attr, 
			function(x) identical(species.info[x[["reactants"]],"ChEBI"],species.info[x[["products"]],"ChEBI"])
	)
	
	spontaneous = sapply( V(reaction.graph)[no.gene]$attr,
			function(x) grepl("spontaneous", x[["name"]], ignore.case=T)
	)
	
	new.reaction.graph = deleteReconnect(reaction.graph, V(reaction.graph)[no.gene][translocation | spontaneous])
	
	if(missing.genes & length(V(new.reaction.graph)[no.gene]>0)){ 
		new.reaction.graph = deleteReconnect(new.reaction.graph, 
				V(new.reaction.graph)[no.gene], reconnect.threshold)
	}		
	return(new.reaction.graph)
}

deleteReconnect <- function(graph, delete.vids, reconnect.threshold=vcount(graph), copy.attr=FALSE){
	if(length(delete.vids)==0){stop("Deleted vertices are not provided")}
	V(graph)$delete = FALSE 
	V(graph)[delete.vids]$delete = TRUE
	
	#A subgraph only including vertices to be deleted and their neighbours.
	graph.sub = subgraph.edges(graph,
			E(graph)[V(graph)[delete] %--% V(graph)[nei(V(graph)[delete])|delete]])
	#Calculate the shortest path length in this network between "neighbours" (retained) vertices. 
	#In this graph, the shortest path will have to go through deleted vertices only.
	#Path lengths indicate the number of deleted vertices lying between 2 retained vertices.
	short.paths = shortest.paths(graph.sub, V(graph.sub)[!delete],V(graph.sub)[!delete], "out")
	new.edges = which(short.paths!=0 & short.paths!=Inf & short.paths< reconnect.threshold+1, arr.ind=T)
	new.edges = rbind(rownames(new.edges),colnames(short.paths)[new.edges[,2]])  #edges to be added
	
	#Get Actual shortest paths for vertices that passed the threshold, to copy edge attributes.
	attr = NULL
	if(copy.attr){
		paths = mapply(function(from, to) get.shortest.paths(rgsub, from, to, mode="out"), new.edges[1,], new.edges[2,])
		attr = lapply(paths, function(x) E(rgsub, path=x)$attr)
	}
	
	new.graph = add.edges(graph, new.edges, attr=attr)
	new.graph = delete.vertices(new.graph, V(graph)[delete])
	remove.vertex.attribute(new.graph, "delete")
	return(new.graph)
	
}


assignEdgeWeights <- function(microarray,gene.graph,y = NULL,corRepeats = 100, same.gene.penalty = -1) {
	microarray = as.matrix(microarray)  #Microarray should be a Dataframe or a matrix, with genes as rownames, and samples as columns.
	
	network.genes <- unique( sub("#(.+)$","", V(gene.graph)$name) ) #break the gene:reaction edge labels, get unique list of genes.
	
	intersect.genes <- intersect(network.genes, rownames(microarray))
	
	cat(sum(!rownames(microarray) %in% intersect.genes), "genes were present in the microarray, but not represented in the network.\n")
	cat(sum(!network.genes %in% intersect.genes), "genes were couldn't be found in microarray.\n")
	
			
	# Get the edgelist, 
	edgelist <- sub("#(.+)$","",get.edgelist(gene.graph) )
	edgelist <- match(edgelist, rownames(microarray))-1
	
	
	# correlation function
	compCor <- function(microarray,edgelist,samegene,weight, corRepeats) {
		all.cors <- .C("corEdgeWeights",
				as.double(t(microarray)),
				as.integer(edgelist),
				as.integer(samegene),
				weight = as.double(weight),
				as.integer(length(edgelist)/2),
				as.integer(ncol(microarray)),
				as.integer(corRepeats),
				NAOK = TRUE)
		return(all.cors$weight)(all.cors$weight)
	}
			
		
	edge.weights <- NULL
	if (is.null(y)) {
		cat("Assigning edge weights.\n")      
		edge.weights <- compCor(microarray,edgelist,E(gene.graph)$samegene,double(length(edgelist)/2), corRepeats)
		edge.weights <- data.frame(edge.weights)
		names(edge.weights) <- "cor"
		
	} else {
		for (yl in 1:nlevels(y)) {  
			data <- microarray[,y == levels(y)[yl]]
			cat("Assigning edge weights for label",levels(y)[yl], "\n") 
			ed <- compCor(microarray,edgelist,E(gene.graph)$samegene,double(length(edgelist)/2), corRepeats)
			if (is.null(edge.weights)) edge.weights <- data.frame(ed)
			else edge.weights <- cbind(edge.weights,data.frame(ed))
		}
		names(edge.weights) <- paste("cor:",levels(y),sep = "")
		
	}
	
	edge.weights[E(gene.graph)$samegene,] <- same.gene.penalty
	E(gene.graph)$edge.weights = as.list(as.data.frame(t(edge.weights)))
	gene.graph$y.labels = if(is.null(y)) "" else levels(y)
	return(gene.graph)
}
