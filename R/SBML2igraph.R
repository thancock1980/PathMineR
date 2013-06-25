# TODO: Add comment
# 
# Author: ahmedmohamed
###############################################################################


SBML2igraph <- function(filename){
	if (!file.exists(filename)) stop("Cannot find file!")
	
	zsbml <- .Call("readsbmlfile", FILENAME = filename)
	
	print("SBML file processed successfully :|")
	print("Processing Gene Association file")
	#################Gene Association File ###################
	galines <- read.delim("ga.txt", header=FALSE, stringsAsFactors=FALSE, skip=1)
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
	
	print("Constructing the iGraph Object :||")
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
		
	reactions <- sapply(V(graph)[type]$name, function(x) reaction.info[x], USE.NAMES=FALSE)
	V(graph)[names(reactions)]$attr <- reactions
	V(graph)[type]$genes = lapply(V(graph)[type], "[[", "genes") 
	
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

makeReactionNetwork <- function(graph){
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
	V(reaction.graph)$size <- 5
	E(reaction.graph)$arrow.size <- 0.2
	E(reaction.graph)$attr <- apply(reaction.edges, 1, 
			function(x) reaction.connections[reaction.connections[,1]==x[1] & reaction.connections[,3]==x[2], 2])
	
	reaction.graph$species.info <- do.call("rbind", lapply(V(graph)[type==FALSE]$attr, "data.frame", stringsAsFactors=FALSE))	
	reaction.graph$gene.list <- graph$gene.list	
	reaction.graph$reaction.coor = cbind(reaction=1:vcount(reaction.graph), 
							met.reaction=match(V(reaction.graph)$name, V(graph)$name) )
#	reaction.graph$reaction.connections <- reaction.connections
	
	return(reaction.graph)
}

makeGeneNetwork.old <- function(reaction.graph){
	genes <- lapply(V(reaction.graph)$attr, "[[","genes")
	names(genes) <- V(reaction.graph)$name
	edges <- get.edgelist(reaction.graph)
	
	#Substitute reaction edges with the corresponding proteins, reconstruct the network
	gene.connections <- do.call("rbind" ,lapply(1:length(edges[,1]), function(x) expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)))
	
	names(gene.connections) <- c("from","to", "idx")
	# edges with the same start and end
	rmidx <- which(gene.connections$from == gene.connections$to)
	if (length(rmidx) > 0) gene.connections <- gene.connections[-rmidx,]
	gene.connections <- unique(gene.connections)
	
	gene.edges <- gene.connections[,-3]
	gene.edges <- unique(gene.edges)
	
	#Make the iGraph object
	gene.graph <- graph.data.frame(gene.edges)
#	V(gene.graph)$id <- gnds
	V(gene.graph)$attr <- reaction.graph$gene.list[V(gene.graph)$name]
	V(gene.graph)$size <- 5
	E(gene.graph)$arrow.size <- 0.2
	
	#Constructing edge attributes
	substrates = E(reaction.graph)[gene.connections$idx]$attr
	names(substrates) = do.call("paste", c(as.data.frame(edges), sep="->") )[gene.connections$idx]
	attr = apply(gene.edges, 1, 
			function(x) substrates[which(gene.connections[,1]==x[1] & gene.connections[,2]==x[2])])
	
	E(gene.graph)$attr <- attr

	
	#Exporting Graph attributes
	gene.graph$species.info <- reaction.graph$species.info
	reaction.info <- V(reaction.graph)$attr
	names(reaction.info) <- V(reaction.graph)$name
	gene.graph$reaction.info <- reaction.info
	
	return(gene.graph)
}

makeGeneNetwork <- function(reaction.graph){
	edges <- get.edgelist(reaction.graph)
	genes <- lapply(V(reaction.graph)$attr, "[[","genes")
	names(genes) <- V(reaction.graph)$name
	genes = mapply(paste, genes, names(genes))
	
	# Convert reaction edges to gene edges.
	gene.connections <- do.call("rbind" ,
			lapply(1:dim(edges)[1], function(x) 
							expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)
	))
	
	#make the graph.
	gene.graph = graph.data.frame(gene.connections[,-3])
	V(gene.graph)$attr <- reaction.graph$gene.list[sub(" \\w+$", "",V(gene.graph)$name)]
	V(gene.graph)$size <- 5
	E(gene.graph)$arrow.size <- 0.2
	
	# get the edge attributes
	E(reaction.graph)$attr[gene.connections$Var3]
	
	# Map vertex indices to the reaction graph (for the layout function). 
	gene.coor <- cbind(gene=1:vcount(gene.graph), 
					reaction=match(sub("^\\w* ", "",V(gene.graph)$name), V(reaction.graph)$name))
	# get graph attributes.
	gene.graph$species.info <- reaction.graph$species.info
	reaction.info <- V(reaction.graph)$attr
	names(reaction.info) <- V(reaction.graph)$name
	gene.graph$reaction.info <- reaction.info
	
	return(gene.graph)
}

makeALL <- function(){
	data("ALL")
	genotype <- which(as.character(ALL$mol.biol) %in% c("NEG", "BCR/ABL"))
	bcell <- grep("^B", ALL$BT)
	ALLsub <- ALL[, intersect(bcell,genotype)]
	
	return(ALLsub)
}

getExpression <- function(ALLsub){
	probe.names <- unique(row.names(exprs(ALLsub)))
	
	mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	uniprot.acc <- getBM( attributes=c("uniprot_swissprot_accession", "affy_hg_u95av2"), 
					filters="affy_hg_u95av2", values=probe.names, mart=mart)
		
	uniprot.acc <- uniprot.acc[uniprot.acc[,1]!="", ]	
}

getProbeInfo <- function(genes, microarray){
	probe.gene.idx <- lapply(genes, function(x) match(uniprot.acc[uniprot.acc[,1]==x, 2], row.names(microarray)))
	names(probe.gene.idx) <- genes
	return(probe.gene.idx)	
}

computeMaxCorrelation <- function(reaction.graph){
	print("1")
	print(proc.time())
	
	load("pathscope.RData")	
	microarray <- exprs(ALLsub)
	
	genes <- lapply(V(reaction.graph)$attr, "[[","genes")
	names(genes) <- V(reaction.graph)$name
	edges <- get.edgelist(reaction.graph)
	
	print("1^")
	print(proc.time())
	
	#Substitute reaction edges with the corresponding proteins, reconstruct the network
	gene.connections <- do.call("rbind" ,lapply(1:length(edges[,1]), function(x) expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)))
	
	
	
	names(gene.connections) <- c("from","to", "idx")
	# edges with the same start and end
	rmidx <- which(gene.connections$from == gene.connections$to)
	if (length(rmidx) > 0) gene.connections <- gene.connections[-rmidx,]
	gene.connections <- unique(gene.connections)
	
	gene.edges <- gene.connections[,-3]
	gene.edges <- unique(gene.edges)
	
	print("1^^")
	print(proc.time())
	
	#mark reverse edges
	forward.edges <-  do.call("paste", c(gene.edges, sep="->") )
	reverse.edges <- do.call("paste", c(gene.edges[,c(2,1)], sep="->") )
	repeated.edges <- match(forward.edges, reverse.edges)
	
	# Mark all edges to be computed, leave out duplicate edges, with their reverse indices
	repeated.edges[-which(1:length(repeated.edges) > repeated.edges)] <- -1
	
	# Get probe ids for genes.
	unique.genes <- unique(unlist(gene.edges))
	probe.gene.idx <- getProbeInfo(unique.genes, microarray)
	
	print("1^^^")
	print(proc.time())
	
	# Mark edges with genes that doesn't have matched probe ids.
	no.probe = names(Filter(function(x) length(x)<1, probe.gene.idx))
	repeated.edges[which(gene.edges[,1]%in% no.probe | gene.edges[,2]%in% no.probe)] <- -2
	
	print("2")
	print(proc.time())
	# Compute Pearson correlation
	cor = lapply(1:length(gene.edges[,1]), 				# For All edges
	function(x) ifelse(repeated.edges[x]==-1,		# Edges marked for computation (excluding duplicate, no expression edges)
			correlateExpression(microarray[ unlist(probe.gene.idx[[ gene.edges[x,1] ]]) ,],  
					microarray[ unlist(probe.gene.idx[[ gene.edges[x,2] ]]),] , 
					multiple.probe.method="mean", bootstrap=100), 
			-2))									# Otherwise, set correlation as -2
	
	# Assign correlations to reverse edges
	cor[ which(repeated.edges>-1) ] = cor[ repeated.edges[repeated.edges>-1] ]
	
	print("3")
	print(proc.time())
	# Assign maximum correlation to Reaction edges
	# Currently indices to reaction edges are stored in "idx" column in "gene connections"
	# Correlations are computed for "gene edges", which is a non-redudant version of "gene connections"
	forward.connections = do.call("paste", c(as.list(gene.connections[, c(1,2)]), sep="->") )
	con2edge = match(forward.connections, forward.edges)		#index table to convert from gene connections to edges.
	
	reaction.edge.cor <- sapply(1:length(edges[,1]),		#be careful, some Rx don't have any gene conn. (cor= -Inf)
		function(x) max(unlist(  cor[ con2edge[ which(gene.connections[,3]==x) ]] )) )
	
	
	E(reaction.graph)$cor <- reaction.edge.cor 
	edge.color <- colorRampPalette(c("light green", "yellow", "orange", "red"))
	E(reaction.graph)$color <- ifelse( E(reaction.graph)$cor < -1, "gray", edge.color(20) )
	
	print("4")
	print(proc.time())
	
	# Convert gene edges to numerical indices
#	gene.edges.num <- cbind(match(gene.edges[,1], names(probe.gene.idx)),
#						match(gene.edges[,2], names(probe.gene.idx)))
#					
#	cors <- .Call("compute_edge_correlation", gene.edges.num, repeated.edges, as.matrix(microarray), probe.gene.idx)
#	
	return(reaction.graph)
}

correlateExpression <- function(probeset1, probeset2, multiple.probe.method="mean", bootstrap = 1){
		
	gene1 = if(is.vector(probeset1)){
		gene1 = probeset1		
	} else{ 
		if(multiple.probe.method=="mean"){
			gene1 = colMeans(probeset1)
		}else{ gene1 = apply(probeset1, 2, "median")}
	}
	
	gene2 = if(is.vector(probeset2)){
		gene2 = probeset2		
	} else{ 
		if(multiple.probe.method=="mean"){
			gene2 = colMeans(probeset2)
		}else{ gene2 = apply(probeset2, 2, "median")}
	}
	w <- length(gene1)	
	s = sample(w, w, replace=TRUE)
	
	c.list = lapply(1:bootstrap, function(x) cor(gene1[s], gene2[s], method="pearson"))
		
	return(median(unlist(c.list)))
}