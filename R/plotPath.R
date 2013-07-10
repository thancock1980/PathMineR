plotPath <- function(graph, path){
	
	#fixing layout start and end nodes
	minx = rep(-1, vcount(graph));	maxx = rep(1, vcount(graph));	miny = rep(-1, vcount(graph));	maxy = rep(1, vcount(graph))
	minx[c(path[1],path[length(path)])] = c(-1,1); maxx[c(path[1],path[length(path)])] = c(-1,1)
	miny[c(path[1],path[length(path)])] = 0; maxy[c(path[1],path[length(path)])] = 0
	
	l = layout.fruchterman.reingold(graph, minx = minx*100, maxx = maxx*100, miny = miny*100, maxy = maxy*100, niter=100000)
	graph.mark = get.edgelist(graph, F)[E(graph, path=path),]
	
	
	reaction.graph = makeReactionNetwork(graph)
	reaction.coor = match(V(reaction.graph)$name,V(graph)$name)
	l2 = l[reaction.coor,]
	reaction.path = as.vector(na.omit(match(V(graph)$name,V(reaction.graph)$name)[path]))
	reaction.mark = get.edgelist(reaction.graph, F)[E(reaction.graph, path=reaction.path),]
	

	
	# make Gene graph
	#ggraph = makeGeneGraph(rgraph)
	edges <- get.edgelist(reaction.graph)
	genes <- lapply(V(reaction.graph)$attr, "[[","genes")
	names(genes) <- V(reaction.graph)$name
	genes = mapply(paste, genes, names(genes))
	gene.connections <- do.call("rbind" ,lapply(1:length(edges[,1]), function(x) expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)))
	gene.graph = graph.data.frame(gene.connections[,1:2])
	V(gene.graph)$size <- 5
	E(gene.graph)$arrow.size <- 0.2
	
	#Fix Coordinates of gene graph
	gene.coor = match(sapply(genes, function(x) x[1]), V(gene.graph)$name)
	minx = rep(-100, vcount(gene.graph));	maxx = rep(100, vcount(gene.graph));	miny = rep(-100, vcount(gene.graph));	maxy = rep(100, vcount(gene.graph))
	minx[gene.coor] = l2[,1]; maxx[gene.coor] = l2[,1]
	miny[gene.coor] = l2[,2]; maxy[gene.coor] = l2[,2]
	l3 = layout.fruchterman.reingold(gene.graph, minx = minx, maxx = maxx, miny = miny, maxy = maxy, niter=100000)
	
	gene.path = gene.coor[reaction.path]
	gene.mark = get.edgelist(gene.graph, F)[E(gene.graph, path=gene.path),]
	
	par(mfrow=c(1,3), mar=c(0,0,5,0)+0.1)
	plot(graph, layout=l, vertex.label="", main="Metbolite-Reaction Network", frame=T,
			mark.groups=split(graph.mark,row(graph.mark)), mark.col=rgb(t(col2rgb("lightgreen")), alpha=100,maxColorValue=255), mark.border=NA, mark.shape=0, mark.expand=10)
	plot(reaction.graph, layout=l2, vertex.label="", main="Reaction-Reaction Network", frame=T,
			mark.groups=split(reaction.mark,row(reaction.mark)), mark.col="lightgreen", mark.border=NA, mark.shape=0, mark.expand=10)
	plot(gene.graph, layout=l3, vertex.label="", main="Gene Network", frame=T,
			mark.groups=split(gene.mark,row(gene.mark)), mark.col="lightgreen", mark.border=NA, mark.shape=0, mark.expand=10)
}