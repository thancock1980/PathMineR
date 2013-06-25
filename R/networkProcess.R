getKEGG <- function(maps = "all",organism = "",kegg.xml.location = "",url = "ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/organisms") {
  mergedKEGG <- NULL
  hasXML <- require(XML)
  if (hasXML == FALSE) stop("To allow parsing of the KGML files please install the XML package.")

  curl <- NULL
  if (kegg.xml.location == "") {
    if (require(RCurl)) {
        cat("\nOpening connection to:",paste(url,"/",organism,"/",sep = ""))
        curl <- getCurlHandle()
    } else stop("To enable downloading of the KGML files please install the RCurl package.")
  }

  # Process the file names so that everything is expected
  if (maps[1] == "all") {
    # get all kegg pathway maps on the ftp site
    if (kegg.xml.location == "") maps <- getKEGGFileNames(organism,url = url)
    else {
      # get all kegg pathway maps in a local directory
      maps <- list.files(kegg.xml.location)
      maps <- grep(maps,pattern = ".xml",value = TRUE)
    }
    # remove the map corresponding to complete metabolism because we use the pathway groups in the annotations
    maps <- maps[-grep("01100",maps)]
    maps <- maps[-grep("01110",maps)]

  } else {
    # ensure the specified maps have a xml extension
    idx <- grep(maps,pattern = ".xml")
    if (length(idx) != length(maps)) {
      if (length(idx) > 0) maps[-idx] <- paste(maps[-idx],".xml",sep = "")
      else maps <- paste(maps,".xml",sep = "")
    }
  }

  cat("\nProcessing KEGG XML files")
  # get the XML files
  for (m in maps) {
    KEGGxml <- NULL

    if (kegg.xml.location == "") {   
      cat("\n  Downloading KEGG XML map:",m)
      xmlTxt <- getURL(paste(url,"/",organism,"/",m,sep = ""),curl = curl)
      KEGGxml <- xmlTreeParse(xmlTxt)
    } else {
      cat("\n  Loading KEGG XML map:",m)
      KEGGxml <- xmlTreeParse(paste(kegg.xml.location,"/",m,sep = ""))
    }

    cat("\n  Converting XML to compound network")
    m.network <- getKEGGInfo(KEGGxml)

    if (!is.null(m.network)) {
      cat("\n  Merging into the overall compound network")
      if (is.null(mergedKEGG)) {
        mergedKEGG <- m.network
      } else {
        # Merge the compound network
        mergedKEGG$edge.list <- rbind(mergedKEGG$edge.list,m.network$edge.list)
        idx <- which(duplicated(mergedKEGG$edge.list))
        if (length(idx) > 0) mergedKEGG$edge.list <- mergedKEGG$edge.list[-idx,]

        # Merge the reaction list
        idx <- which(names(m.network$reaction.info) %in% names(mergedKEGG$reaction.info))
        if (length(idx) > 0) {
          for (r in idx) {
            g1 <- mergedKEGG$reaction.info[[ names(m.network$reaction.info)[r] ]]
            g2 <- m.network$reaction.info[[r]]
            g <- unique(c(g1,g2)) # if a reaction occurs in both network, take the gene set union
            mergedKEGG$reaction.info[[ names(m.network$reaction.info)[r] ]] <- g
          }
          mergedKEGG$reaction.info <- c(mergedKEGG$reaction.info,m.network$reaction.info[-idx])
        } else mergedKEGG$reaction.info <- c(mergedKEGG$reaction.info,m.network$reaction.info)
        mergedKEGG$details <- rbind(mergedKEGG$details,m.network$details)
      }
    }
  }

   return(mergedKEGG)
}

getKEGGFileNames <- function(organism,url) {
  kfiles <- c()

  organism.url <- paste(url,"/",organism,"/",sep = "")
  cat("\nGetting KEGG maps from",organism.url)
  
  kegginfo <- getURL(organism.url)

  kfiles <- strsplit(strsplit(kegginfo,"\n")[[1]]," ")
  kfiles <- unlist(sapply(kfiles,grep,pattern = ".xml",value = TRUE))

  return(kfiles)
}

# Takes the location and code of a KEGG XML file and finds a reaction->compound edge list and extracts
# which genes belong to each reaction
getKEGGInfo <- function(KEGG.XML) {
	pathway <- KEGG.XML$doc$children$pathway
	pathtitle <- xmlGetAttr(pathway,"title")
    pathcode <- strsplit(xmlGetAttr(pathway,"name"),":")[[1]][2]
    if (is.null(pathtitle)) pathtitle <- xmlGetAttr(pathway,"name")

	reaction.info <- list()
	edges <- data.frame(from = c(),to = c(),label = c(),pathway = c(),stringsAsFactors = FALSE)

	r.indx <- which(names(pathway) %in% "reaction")
	e.indx <- which(names(pathway) %in% "entry")
	if (length(r.indx) == 0) {
        cat("\n    !!!WARNING: No reactions found in pathway map!!!")
		return(NULL)
    }

	for (reaction in r.indx) {
        all.r.names <- strsplit(xmlGetAttr(pathway[[reaction]],"name")," ")[[1]]
        all.r.names <- all.r.names[all.r.names != ""]
        # One compound-compound transfer can occur by more than one reaction
        for (r.name in all.r.names) {
            disp.r.name <- strsplit(r.name,split = ":")[[1]][2]

            # First look to see if the new reaction has any valid genes
            reaction.has.genes <- FALSE
            for (entry in e.indx) {
                e.reaction <- xmlGetAttr(pathway[[entry]],"reaction")      
                if (!is.null(e.reaction)) { 
                    e.reaction <- strsplit(e.reaction,split = " ")[[1]]
                    e.reaction <- e.reaction[e.reaction != ""]
                    if (r.name %in% e.reaction) { 
                        reaction.genes <- strsplit(xmlGetAttr(pathway[[entry]],"name"),split = " ")[[1]]
                        reaction.genes <- reaction.genes[reaction.genes != ""] # remove whitespace genes
                        reaction.genes <- strsplit(toupper(reaction.genes),":")
                        
                        # exclude all KEGG ortholog elements as these wont be found in a microarray
                        ko.genes <- which(sapply(reaction.genes,"[[",1) == "KO") 
                        if (length(ko.genes) > 0) reaction.genes <- reaction.genes[-ko.genes]
                        reaction.genes <- sapply(reaction.genes,"[[",2)

                        if (length(reaction.genes) > 0) {
                            reaction.has.genes <- TRUE
                            if (disp.r.name %in% names(reaction.info)) {
                                reaction.info[[disp.r.name]] <-    
                                    c(reaction.info[[disp.r.name]],reaction.genes)
                            } else {
                                reaction.info <- c(reaction.info,list(reaction.genes))
                                names(reaction.info)[length(reaction.info)] <- disp.r.name
                            }
                        }
                    }
                }
            }
            
            # if the reaction has genes then add it to the network
            if (reaction.has.genes) {
                # add edge to compound-compound network labelled by reaction
                r.type <- xmlGetAttr(pathway[[reaction]],"type")
                
                r.from <- which(names(pathway[[reaction]]) %in% "substrate")
                r.to <- which(names(pathway[[reaction]]) %in% "product")
                for (sub in r.from) {
                    r.from.name <- strsplit(xmlGetAttr(pathway[[reaction]][[sub]],"name"),split = ":")[[1]][2]
                    for (prod in r.to) {
                        r.to.name <- strsplit(xmlGetAttr(pathway[[reaction]][[prod]],"name"),split = ":")[[1]][2]
      
                        edges <- rbind(edges,
                            data.frame(from = r.from.name,to = r.to.name,label = disp.r.name,pathway = pathtitle,
								   stringsAsFactors = FALSE))

                        if (r.type == "reversible") 
                            edges <- rbind(edges,
                                data.frame(from = r.to.name,to = r.from.name,label = disp.r.name,pathway = pathtitle,
									  stringsAsFactors = FALSE))
                    }
                } 
            }
		}
	}

	return(list(edge.list = edges,reaction.info = reaction.info,
                details = data.frame(code=pathcode,title = pathtitle,stringsAsFactors = FALSE)))
}

gannotate <- function(gene,from,to,reaction,pathway,ann = "GRFT") {
  if (ann == "GRPFT") gl <- paste(gene,reaction,pathway,from,to,sep=":")
  if (ann == "GRFT") gl <- paste(gene,reaction,from,to,sep=":")
  if (ann == "GRFT") gl <- paste(gene,pathway,from,to,sep=":")
  if (ann == "GFT") gl <- paste(gene,from,to,sep=":")
  if (ann == "GR") gl <- paste(gene,reaction,sep=":")
  return(gl) 
}

removeUnobservedGenes <- function(kobj,microarray=NULL) {
  compoundNetwork <- kobj[[1]]
  reactionList <- kobj[[2]]
 
  micgenes <- unique(names(microarray))
  rmreact <- c()
  for (i in 1:length(reactionList)) {
    if (length(reactionList[[i]]) == 0) {
      rmreact <- c(rmreact,i) #empty reactions
    } else if (!is.null(micgenes)) {       
      rmidx <- which((reactionList[[i]] %in% micgenes) == FALSE) # Unobserved genes
      if (length(rmidx) > 0) {
        reactionList[[i]] <- reactionList[[i]][-rmidx]
        if (length(reactionList[[i]]) == 0) rmreact <- c(rmreact,i)
      }
    }
  }
  if (length(rmreact) > 0) {
    compoundNetwork <- compoundNetwork[-which(compoundNetwork$label %in% names(reactionList)[rmreact]),]
    reactionList <- reactionList[-rmreact]
  }

  kobj$edge.list <- compoundNetwork
  kobj$reaction.info <- reactionList

  return(kobj)
}

# take the output of kegg.info and make it a gene network for ichigaku's model.
computeGeneNetwork <- function(kegg.obj,ann = "GFTR") {
  cat("\nConverting compound/reaction network to gene/compound network")

  compoundNetwork <- kegg.obj[[1]]
  reactionList <- kegg.obj[[2]]

  # Gene label = GENE_NAME:CFROM:CTO;REACTION:PATHWAY
  # this unfortunately is the most concise descriptor for every possible node in KEGG such that
  # broken pathways done occur.
  pstfunc <- function(i,compoundNetwork,reactionList) {
    gl <- reactionList[[compoundNetwork[i,]$label]]
    gl <- gannotate(gl,compoundNetwork[i,]$from,compoundNetwork[i,]$to,compoundNetwork[i,]$label,compoundNetwork[i,]$pathway,ann)
    return(gl)
  }
  genelist <- lapply(1:nrow(compoundNetwork),pstfunc,compoundNetwork,reactionList)
  genelist <- unique(unlist(genelist))

  # Convert a compound network into a gene network *brackets indicate the edge label
  # Compound network: C1 --> (g1) --> C2 --> (g2) --> C3
  # Gene network: g1 --> (C1:C2:C3) --> g2
  geneNetwork <- NULL
  for (i in 1:nrow(compoundNetwork)) {
    # process the first set start compound -> middle
    cfrom <- as.character(compoundNetwork[i,]$from)
    rfrom <- as.character(compoundNetwork[i,]$label) # reactions joining start and middle compounds
    gl <- gannotate(reactionList[[rfrom]],cfrom,compoundNetwork[i,]$to,rfrom,compoundNetwork[i,]$pathway,ann)
    gfrom <- which(genelist %in% gl)  # index of genes joining the start and middle compounds
    
    # process the next step, middle compound -> next set of compounds
    middle <- as.character(compoundNetwork[i,]$to) # middle compound - also the edge label
    path <- as.character(compoundNetwork[i,]$pathway)
    edgeIdx <- which(compoundNetwork$from == middle) # indices reactions steming from the middle compound
    if (length(edgeIdx) > 0) { # if the middle compound is not a terminating compound
      zcomps <- compoundNetwork[edgeIdx,]
      idx <- which(zcomps$to == cfrom) # cannot have a self loop from cfrom -> middle -> cfrom
      if (length(idx) > 0) zcomps <- zcomps[-idx,] 
      
      if (nrow(zcomps) > 0) {
               edges <- NULL
               for (j in 1:nrow(zcomps)) {
                    rto <- as.character(zcomps[j,]$label) # reactions joining middle and end compounds
                    path <- zcomps[j,]$pathway
                    cend <- zcomps[j,]$to
                    gto <- length(genelist)
                    gl <- gannotate(reactionList[[rto]],zcomps[j,]$from,
                                      zcomps[j,]$to,zcomps[j,]$label,zcomps[j,]$pathway,ann)
                    gto <- which(genelist %in% gl) # get the gene indices of gto

                    # make the edges and append to the geneNetwork
                    # interaction always throws a warning, all this does is
                    # get all combinations of gfrom and gto which may be of different sizes
                    eds <- suppressWarnings(levels(interaction(gfrom,gto))) # create all edge combinations of gfrom -> gto
                    eds <- lapply(strsplit(eds,".",fixed = TRUE),as.numeric)
                    eds <- data.frame(t(do.call("cbind",eds)))
                    idx <- which(eds[,1] == eds[,2]) #remove edges from one gene to itself
                    if (length(idx) > 0) eds <- eds[-idx,]

                    if (nrow(eds) > 0) {
                        # Each edge is labelled by C_from:C_middle:C_end and pathway (C = compound)
                        # Each gene is labelled by Gene:Reaction
                        # So the edge information is complete
                        label <- paste(c(cfrom,middle,cend),collapse = ":")
                        eds <- data.frame(cbind(eds,label,path))
                        names(eds) <- c("from","to","label","pathway")
                        if (is.null(edges)) edges <- eds
                        else edges <- rbind(edges,eds)
                    }
               }

               # if there are any edges left then add them to the network
               if (!is.null(edges)) {                  
                    if (is.null(geneNetwork)) geneNetwork <- edges
                    else geneNetwork <- rbind(geneNetwork,edges)
               }
      }
    }
  }

  idx <- which(duplicated(geneNetwork))
  if (length(idx) > 0) geneNetwork <- geneNetwork[-idx,]
  geneNetwork <- with(geneNetwork,geneNetwork[order(from,to),])
  dimnames(geneNetwork)[[1]] <- 1:nrow(geneNetwork)
  geneNetwork <- list(nodes = genelist,edges = geneNetwork)

  geneNetwork$edges$to <- as.integer(geneNetwork$edges$to)
  geneNetwork$edges$from <- as.integer(geneNetwork$edges$from)
  geneNetwork$edges$label <- as.character(geneNetwork$edges$label)
  geneNetwork$edges$pathway <- as.character(geneNetwork$edges$pathway)

  return(geneNetwork)
}

defineNeighborhood <- function(kobj,kgn,edge.weights,seed,range,ann = "GRFT",mode = "out") {
    nobj <- kgn
    nobj$nodes <- c("s",kgn$nodes,"t")
    nobj$edges[,1:2] <- nobj$edges[,1:2] + 1

    # new start compounds
    sout <- NULL
    sedges <- kobj$edge[kobj$edge$from %in% seed,]
    if (nrow(sedges) == 0) {
        cat("\nSeed Compound",seed,"No Start Edges")
        return(NULL)
    }
    sedges <- as.character(apply(sedges,1,FUN = paste,collapse = ":"))
    for (i in 1:length(sedges)) {
        s <- strsplit(sedges[i],":")[[1]]
        gnms <- gannotate(kobj$reaction.info[[s[3]]],s[1],s[2],s[3],s[4],ann = ann)
        gidx <- which(nobj$nodes %in% gnms)
        gl <- length(gidx)
        se <- data.frame(from = rep(1,gl),
                         to = gidx,
                         label = rep(paste(c("s",s[1],s[2]),collapse = ":"),gl),
                         pathway=rep(s[4],gl),stringsAsFactors=FALSE)
        if (is.null(sout)) sout <- se
        else sout <- rbind(sout,se)
    }
    sout <- unique(sout)   

    # get the neighborhood
    e <- kobj$edge[,1:2]
    vlist <- unique(c(e$from,e$to))
    em <- data.frame(from = match(e$from,vlist),to = match(e$to,vlist))
    ge <- graph.edgelist(el = as.matrix(em),directed = TRUE)
    s <- which(vlist == seed)
    nidx <- neighborhood(ge,order = range,nodes = s,mode = mode)[[1]]
    if (length(nidx) == 0) {
        cat("\nSeed Compound",seed,"has no neighborhood")
        return(NULL)
    }
    end <- sort(vlist[nidx[nidx != s]])

    # new end compounds
    tout <- NULL
    tedges <- kobj$edge[kobj$edge$to %in% end,]
    if (nrow(tedges) == 0) {
        cat("\nSeed Compound",seed,"No End Edges")
        return(NULL)
    }
    tedges <- as.character(apply(tedges,1,FUN = paste,collapse = ":"))
    for (i in 1:length(tedges)) {
        t <- strsplit(tedges[i],":")[[1]]
        gnms <- gannotate(kobj$reaction.info[[t[3]]],t[1],t[2],t[3],t[4],ann = ann)
        gidx <- which(nobj$nodes %in% gnms)
        gl <- length(gidx)
        te <- data.frame(from = gidx,
                         to = rep(length(nobj$nodes),gl),
                         label = rep(paste(c(t[1],t[2],"t"),collapse = ":"),gl),
                         pathway="end",stringsAsFactors=FALSE)
        if (is.null(tout)) tout <- te
        else tout <- rbind(tout,te)
    }
    tout <- unique(tout)

    nobj$edges <- rbind(sout,nobj$edges,tout)
    zidx <- order(nobj$edges[,1],nobj$edges[,2])
    nobj$edges <- nobj$edges[zidx,]
    rownames(nobj$edges) <- as.character(1:nrow(nobj$edges))

    sw <- data.frame(matrix(1,nrow = nrow(sout),ncol = ncol(edge.weights)));names(sw) <- names(edge.weights)
    tw <- data.frame(matrix(1,nrow = nrow(tout),ncol = ncol(edge.weights)));names(tw) <- names(edge.weights)
    weights <- data.frame(rbind(sw,edge.weights,tw)[zidx,])
    names(weights) <- names(edge.weights)    

    nobj$edges$to <- as.integer(nobj$edges$to)
    nobj$edges$from <- as.integer(nobj$edges$from)
    nobj$edges$label <- as.character(nobj$edges$label)
    nobj$edges$pathway <- as.character(nobj$edges$pathway)

    return(list(geneNetwork = nobj,edge.weights = weights,seed = seed,neighborhood = end))
}

assignEdgeWeights <- function(microarray,geneNetwork,y = NULL,corRepeats = 100,sameGenePenalty = -1) {
   # correlation function
   compCor <- function(rowSamples,all.genes,from,to,weight,samegene) {
        ag <- all.genes[rowSamples,]
        all.cors <- .C("SCOPE_corEdgeWeights",
            as.double(as.matrix(ag)),
            as.integer(from),
            as.integer(to),
            as.integer(samegene),
            weight = as.double(weight),
            as.integer(length(from)),
            as.integer(nrow(ag)),NAOK = TRUE)
        return(all.cors$weight)(all.cors$weight)
    }  

   gn.r <- geneNetwork$nodes
   gn <- sapply(strsplit(gn.r,":"),"[[",1) #break the gene:reaction edge labels

   Rmicroarray <- data.frame(lapply(microarray,rank),stringsAsFactors=FALSE)
   names(Rmicroarray) <- names(microarray)

   stvars <- data.frame(matrix(0,nrow = nrow(Rmicroarray),ncol = 2))
   names(stvars) <- c("s","t")
   p.microarray <- cbind(microarray,stvars)
   all.genes <- p.microarray[gn] # extract all relevant genes
   names(all.genes) <- gn.r # relabel by the gene:reaction labels

   # Compute all gene-gene correlations
   from <- geneNetwork$edges[,1] - 1
   to <- geneNetwork$edges[,2] - 1 

   # find all same gene edges
   gt <- sapply(strsplit(geneNetwork$nodes[geneNetwork$edges$to],":"),"[[",1)
   gf <- sapply(strsplit(geneNetwork$nodes[geneNetwork$edges$from],":"),"[[",1)
   sgid <- which(gf == gt)
   samegene <- integer(length(gt)); samegene[sgid] <- 1
    
   edge.weights <- NULL
   if (is.null(y)) {
        cat("\nAssigning edge weights")      
        if (corRepeats > 1) {
            cat("\nTaking the median of",corRepeats,"Bootstrapped correlations - please wait a bit...")
            bsize <- ifelse(nrow(all.genes) < 10,10,nrow(all.genes))
            w.list <- lapply(1:corRepeats,function(i) sample(1:nrow(all.genes),bsize,replace = TRUE))
            c.list <- lapply(w.list,compCor,all.genes,from,to,double(length(from)),samegene)
            c.list <- do.call(cbind,c.list)
            edge.weights <- apply(c.list,1,median,na.rm = TRUE) 
        } else edge.weights <- compCor(1:nrow(all.genes),all.genes,from,to,double(length(from)),samegene)
        edge.weights <- data.frame(edge.weights)
        names(edge.weights) <- "cor"
   } else {
        for (yl in 1:nlevels(y)) {  
            data <- all.genes[y == levels(y)[yl],]
            cat("\nAssigning edge weights for label",levels(y)[yl]) 
            ed <- c()
            if (corRepeats > 1) {
                cat("\nTaking the median of",corRepeats,"Bootstrapped correlations")
                bsize <- ifelse(nrow(data) < 10,10,nrow(all.genes))
                w.list <- lapply(1:corRepeats,function(i) sample(1:nrow(data),bsize,replace = TRUE))
                c.list <- lapply(w.list,compCor,data,from,to,double(length(from)),samegene)
                c.list <- do.call(cbind,c.list)
                ed <- apply(c.list,1,median,na.rm = TRUE) 
            } else ed <- compCor(1:nrow(data),data,from,to,double(length(from)),samegene)
            if (is.null(edge.weights)) edge.weights <- data.frame(ed)
            else edge.weights <- cbind(edge.weights,data.frame(ed))
        }
        names(edge.weights) <- paste("cor:",levels(y),sep = "")
   }

   if (length(sgid) > 0) {
        for (i in 1:ncol(edge.weights)) edge.weights[,i][sgid] <- sameGenePenalty
   }

   return(edge.weights)
}


edgeSums <- function(microarray,geneNetwork,sgp = 0) {
   gn.r <- geneNetwork$nodes
   gn <- sapply(strsplit(gn.r,":"),"[[",1) #break the gene:reaction edge labels

   not.included.genes <-  which((gn %in% names(microarray)) == FALSE)
   not.included.genes <- not.included.genes[-which(not.included.genes %in% c(1,length(gn)))]
   if (length(not.included.genes) > 0) {
        cat("\nWarning: Could not find the following genes in the microarray.\nRemoving those edges and storing the new network in the geneNetwork element of the returned object.")
        cat("\n",gn[not.included.genes])

        midx <- c(which(geneNetwork$edges$from %in% not.included.genes),which(geneNetwork$edges$to %in% not.included.genes))
        gn <- gn[-not.included.genes]
        gn.r <- gn.r[-not.included.genes]
        
        geneNetwork$edges <- geneNetwork$edges[-midx,]
        geneNetwork$edges$to <- match(geneNetwork$nodes[geneNetwork$edges$to],gn.r)
        geneNetwork$edges$from <- match(geneNetwork$nodes[geneNetwork$edges$from],gn.r)
        geneNetwork$nodes <- gn.r
   }

   stvars <- data.frame(matrix(0,nrow = nrow(microarray),ncol = 2))
   names(stvars) <- c("s","t")
   p.microarray <- cbind(microarray,stvars)
   all.genes <- p.microarray[gn] # extract all relevant genes
   names(all.genes) <- gn.r # relabel by the gene:reaction labels

   stvars <- data.frame(matrix(0,nrow = nrow(microarray),ncol = 2))
   names(stvars) <- c("s","t")
   p.microarray <- cbind(microarray,stvars)
   all.genes <- p.microarray[gn] # extract all relevant genes
   names(all.genes) <- gn.r # relabel by the gene:reaction labels

   # find all same gene edges
   gt <- sapply(strsplit(geneNetwork$nodes[geneNetwork$edges$to],":"),"[[",1)
   gf <- sapply(strsplit(geneNetwork$nodes[geneNetwork$edges$from],":"),"[[",1)
   diffgene <- gt != gf

   edge.weights <- data.frame(cor = double(nrow(geneNetwork$edges)))
   for (i in 1:nrow(geneNetwork$edges)) {
        if (diffgene) 
            edge.weights[,1][i] <- sum(log(microarray[,geneNetwork$edges[i,]$to] * microarray[,geneNetwork$edges[i,]$from]))
   }

   for (i in 1:ncol(edge.weights)) edge.weights[,i][which(!diffgene)] <- sgp

   return(list(edge.weights = edge.weights,geneNetwork = geneNetwork))
}
