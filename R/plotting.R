plotScopeMST <- function(scopes,mincoverset,title.info = NULL,loc = NULL,num.hubs = 5) {
  getCompoundScopeNetwork <- function(file) {
      load(file)
      zp <- unlist(lapply(lapply(p$paths,"[[","compounds"),strsplit,":"),recursive = FALSE)
      zpl <- sapply(zp,length)
      zp <- data.frame(do.call("rbind",zp[-which(zpl == 1)]),stringsAsFactors = FALSE)
      zl <- as.factor(paste(zp[,1],zp[,2],sep = ":"))
      zw <- summary(zl,maxsum=nlevels(zl))
      zp <- data.frame(from = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",1), 
                       to = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",2),
                       weights = as.numeric(zw),stringsAsFactors = FALSE)
      ge <- graph.edgelist(as.matrix(zp[,1:2]),directed = TRUE)
      E(ge)$weights <- zp[,3]
      return(ge)
  }

  rescale <- function(x,newspan = c(-1,1)) {
      oldrange <- max(x) - min(x); newrange <- max(newspan) - min(newspan)
      scale <- oldrange/newrange; shift <- min(x) - (scale*-1)
      return((x - shift)/scale)
  }  

  # compute the similarity between scopes, and find hubs in each cover set compounds mst
  scd <- NULL
  cgraphlist <- NULL
  for (i in 1:length(mincoverset)) {
    iidx <- mincoverset[i]
    j <- i + 1
    while(j>i & j <= length(mincoverset)) {
      jidx <- mincoverset[j]
      z <- data.frame(from = names(scopes)[iidx],
                      to = names(scopes)[jidx],
                      weights = length(intersect(scopes[[iidx]],scopes[[jidx]])) )
      if (is.null(scd)) scd <- z
      else scd <- rbind(scd,z)
      j <- j + 1
    }
    ge <- getCompoundScopeNetwork(paste(loc,names(scopes)[iidx],".RData",sep = ""))
    w <- E(ge)$weights; w <- -log(rank(w)/length(w))
    mst <- minimum.spanning.tree(ge,weights = w)
    if (vcount(mst) > 25) {
          ng <- sort(degree(mst),decreasing = TRUE,index.return = TRUE)          
          cgraphlist[[i]] <- V(mst)$name[ng$ix[ng$x > 9]]
      } else {
        cgraphlist[[i]] <- V(mst)$name[-which(V(mst)$name == "s")] 
      }
      names(cgraphlist)[i] <- names(scopes)[iidx]
    }

  #set graph up
  mst <- graph.edgelist(as.matrix(scd[,1:2]),directed = FALSE)
  E(mst)$weights <- scd$weights
  V(mst)$ctype <- rep("compound",vcount(mst))

  # add intermediate hub reactions
  newnodes <- unique(as.character(unlist(cgraphlist)))
  inodes <- as.character(unlist(cgraphlist))
  mst <- add.vertices(mst,nv = length(newnodes),attr = list(name = newnodes,ctype = rep("hub",length(newnodes)) ))
  nreps <- as.character(unlist(mapply(rep,names(cgraphlist),sapply(cgraphlist,length))))
  ed <- cbind(nreps,inodes)
  ed <- as.character(t(unlist(ed)))
  ed <- match(ed,V(mst)$name)-1
  mst <- add.edges(mst,edges = ed,attr = list(weights = rep(10,length(inodes))))

  #node properties
  vcol <- rep("green",vcount(mst))
  vcol[V(mst)$ctype == "reaction"] <- "red"
  vcol[V(mst)$ctype == "hub"] <- "orange"
  vname <- rep("",vcount(mst))
  vname[V(mst)$ctype %in% "compound"] <- V(mst)$name[V(mst)$ctype %in% "compound"]
  vname[V(mst)$ctype %in% "hub"] <- V(mst)$name[V(mst)$ctype %in% "hub"]
  vshape <- "circle"
  vlab.col <- ifelse(V(mst)$ctype %in% "compound","blue","black")
  vsize <- rep(2,vcount(mst))
  ssize <- sapply(scopes[mincoverset],length)
  vsize[V(mst)$ctype %in% "compound"] <- ssize/max(ssize)*5 + 3
  hsize <- degree(mst)[V(mst)$ctype %in% "hub"]
  vsize[V(mst)$ctype %in% "hub"] <- hsize/max(hsize)*3 + 2
  
  #edge properties
  ecol <- rep("green",ecount(mst))
  ed <- suppressWarnings(get.edgelist(mst))
  ecol[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"hub"] <- "orange"
  ecol[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"reaction"] <- "red"
  ewidth <- rep(2,ecount(mst))
  ewidth[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"compound"] <- 6
  ewidth[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"hub"] <- 3

  E(mst)$col <- ecol; V(mst)$col = vcol; V(mst)$shape <- vshape;V(mst)$vlab.col <- vlab.col
  V(mst)$size <- vsize; V(mst)$vname <- vname
  mst <- delete.edges(mst,E(mst)[ecol == "green" | ecol == "red"])
  mst <- delete.vertices(mst,V(mst)[vcol == "red"])
  mst <- delete.vertices(mst,V(mst)[degree(mst) == 0])
  
  #plot
  zlayout <- layout.fruchterman.reingold(mst,niter = 10000,circular = TRUE)
  plot(mst,layout = zlayout,
       vertex.label = V(mst)$vname,vertex.size = V(mst)$vsize,vertex.color = V(mst)$col,vertex.shape = V(mst)$vshape,
       vertex.label.color = V(mst)$vlab.col,vertex.frame.color = V(mst)$col,
       edge.color = E(mst)$col)
  
  mtext(title.info$ntitle,side = 3,outer = TRUE,font = 2,
        cex = title.info$title.cex,line = title.info$title.line)  

  return(mst)

}

plotCompoundMST <- function(cscope,text.cex = 0.3,title.info = NULL) {
  rescale <- function(x,newspan = c(-1,1)) {
    oldrange <- max(x) - min(x); newrange <- max(newspan) - min(newspan)
    scale <- oldrange/newrange; shift <- min(x) - (scale*-1)
    return((x - shift)/scale)
  }  

  zp <- unlist(lapply(lapply(cscope$paths,"[[","compounds"),strsplit,":"),recursive = FALSE)
  zpl <- sapply(zp,length)
  zp <- data.frame(do.call("rbind",zp[-which(zpl == 1)]),stringsAsFactors = FALSE)
  zl <- as.factor(paste(zp[,1],zp[,2],sep = ":"))
  zw <- summary(zl,maxsum=nlevels(zl))
  zp <- data.frame(from = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",1), 
                   to = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",2),
                   weights = as.numeric(zw),stringsAsFactors = FALSE)

  ge <- graph.edgelist(as.matrix(zp[1:2]),directed = TRUE)
  E(ge)$weights <- zp$weights

  # make minimum spanning tree
  w <- E(ge)$weights; w <- -log(rank(w)/length(w))
  mst <- minimum.spanning.tree(ge,weights = w)

  ##PLOT!!!
  sidx <- match("s",V(mst)$name)
  w <- E(mst)$weights; sw <- -log(rank(w)/length(w))  
  if (vcount(mst) > 25) {  
    zlayout <- layout.reingold.tilford(mst,niter = 3000,root = sidx-1,circular = TRUE)
    zlayout[,1] <- rescale(zlayout[,1]);zlayout[,2] <- rescale(zlayout[,2])
    zlayout <- layout.graphopt(mst,start = zlayout,niter = 3000,mass = 15,spring.length = 25)
    zlayout[,1] <- rescale(zlayout[,1]);zlayout[,2] <- rescale(zlayout[,2])
  } else {
    zlayout <- rescale(layout.fruchterman.reingold(mst))
  }

  vcol <- rep("red",vcount(mst))
  vcol[sidx] <- "green"
  vname <- V(mst)$name 
  vshape <- "circle"
  vlab.col <- rep("black",vcount(mst))
  ssize <- degree(ge)
  vsize <- ssize/max(ssize)*5 
  vsize[sidx] <- 3
  plot(mst,layout = zlayout,rescale = FALSE,
       vertex.label = vname,vertex.size = vsize,vertex.color = vcol,
       vertex.shape = vshape, vertex.label.color = vlab.col,vertex.frame.color = vcol,
       edge.color = "black",edge.width = 2,edge.arrow.size = 0.5)

  mtext(title.info$ntitle,side = 3,outer = TRUE,font = 2,
        cex = title.info$title.cex,line = title.info$title.line)

  return(ge)
}
           
plotCompoundScopeTree <- function(cscope,text.cex = 0.3) {
  zp <- unlist(lapply(lapply(cscope$paths,"[[","compounds"),strsplit,":"),recursive = FALSE)
  zpl <- sapply(zp,length)
  zp <- data.frame(do.call("rbind",zp[-which(zpl == 1)]),stringsAsFactors = FALSE)
  zl <- as.factor(paste(zp[,1],zp[,2],sep = ":"))
  zw <- summary(zl,maxsum=nlevels(zl))
  zp <- data.frame(from = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",1), 
                   to = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",2),
                   weights = as.numeric(zw),stringsAsFactors = FALSE)

  zg <- graph.edgelist(as.matrix(zp[1:2]),directed = TRUE)
  nv <- vcount(zg)

  rescale <- function(x,newspan = c(-1,1)) {
    oldrange <- max(x) - min(x); newrange <- max(newspan) - min(newspan)
    scale <- oldrange/newrange; shift <- min(x) - (scale*-1)
    return((x - shift)/scale)
  }

  ewidth <- zp$weights
  for (vn in V(zg)$name) {
    idx <- which(zp$from == vn)
    if (length(idx) > 0) ewidth[idx] <- ewidth[idx]/sum(ewidth[idx])
  }

  par(mar = c(1,1,1,1))
  snode <- which(V(zg)$name == "s")
  zlayout <- layout.reingold.tilford(zg,circularLogical = FALSE,root = snode-1,weights = ewidth)
  zlayout[,1] <- rescale(zlayout[,1]);zlayout[,2] <- rescale(zlayout[,2])
  vcols <- rep("red",nv);vcols[snode] <- "green"
  vsize <- as.numeric(cut(degree(zg),10))/5; vsize[snode] <- 2
  elabel <- zp$weights
  plot(zg,margin = c(0,0,0,0),rescale = FALSE,ylim = c(-1,1),xlim = c(-1,1),asp = h/w,
        layout = zlayout,
        vertex.size = vsize, vertex.color = vcols,
        vertex.label.family = postscriptFonts("Helvetica"),vertex.label = rep("",nv), 
        edge.arrow.size = 0.3,edge.width = ewidth+0.1)
  vlabels <- V(zg)$name
  vlabels[snode] <- "Start"
  text(x = zlayout[,1][snode],y = zlayout[,2][snode],vlabels[snode], 
        cex = text.cex,col = "blue",font = 2)
  text(x = zlayout[,1][-snode],y = zlayout[,2][-snode],vlabels[-snode],srt = 90, 
        cex = text.cex,col = "blue",font = 2,adj = 0)
  return(zp)
}

bipartiteLayout <- function(el,iter) {
  rescale <- function(x,newspan = c(-1,1)) {
    if (length(x) == 1) return(0)
    oldrange <- max(x) - min(x); newrange <- max(newspan) - min(newspan)
    scale <- oldrange/newrange; shift <- min(x) - (scale*-1)
    return((x - shift)/scale)
  }    

  ## convert to numeric
  nel <- matrix(NA,nrow = nrow(el),ncol = ncol(el))
  nel[,1] <- as.numeric(as.factor(el$from))
  nel[,2] <- as.numeric(as.factor(el$to))
  n1.map <- data.frame(label = levels(as.factor(el$from)),position = 1:nlevels(as.factor(el$from)))
  n2.map <- data.frame(label = levels(as.factor(el$to)),position = 1:nlevels(as.factor(el$to)))

  svalues <- sort(c(nrow(n1.map),nrow(n2.map)))

  #break it up (sort of) 
  smallerpositions <- seq(1,2*svalues[2],2*svalues[2]/svalues[1])

  #which is smaller
  sidx <- 2
  if (nrow(n1.map) < nrow(n2.map)) sidx <-1

  for (iter in 1:iter) {
    nel.old <- nel
    coladjust <- iter %% 2

    if (coladjust == 0) {
      # find the median position of each vertex
      nodemedians <- aggregate(nel[,1],list(as.factor(el[,2])),median)
      n2.map$label <- nodemedians[,1]
      n2.map$position <- rank(nodemedians[,2],ties = "first")

      if (sidx == 2) n2.map$position <- smallerpositions[n2.map$position]

      # overwrite positions of each node
      for (i in 1:nrow(n2.map)) nel[,2][el$to == n2.map$label[i]] <- n2.map$position[i]
    } else {
      # find the median position of each vertex
      nodemedians <- aggregate(nel[,2],list(as.factor(el[,1])),median)
      n1.map$label <- nodemedians[,1]
      n1.map$position <- rank(nodemedians[,2],ties = "first")


      if (sidx == 1) n1.map$position <- smallerpositions[n1.map$position]

      # overwrite positions of each node
      for (i in 1:nrow(n1.map)) nel[,1][el$from == n1.map$label[i]] <- n1.map$position[i]
    }  
  }

  n1.map$position <- rescale(n1.map$position)
  n2.map$position <- rescale(n2.map$position)
  for (i in 1:nrow(n2.map)) nel[,2][el$to == n2.map$label[i]] <- n2.map$position[i]
  for (i in 1:nrow(n1.map)) nel[,1][el$from == n1.map$label[i]] <- n1.map$position[i]
  return(list(n1.pos = n1.map,n2.pos = n2.map,edge.map = nel))
}

plotBipartiteScope <- function(scopes,mincoverset,title.info = NULL,loc = NULL,num.hubs = 5,
                               sbml.info = NULL) {
  getCompoundScopeNetwork <- function(file) {
      load(file)
      zp <- unlist(lapply(lapply(p$paths,"[[","compounds"),strsplit,":"),recursive = FALSE)
      zpl <- sapply(zp,length)
      zp <- data.frame(do.call("rbind",zp[-which(zpl == 1)]),stringsAsFactors = FALSE)
      zl <- as.factor(paste(zp[,1],zp[,2],sep = ":"))
      zw <- summary(zl,maxsum=nlevels(zl))
      zp <- data.frame(from = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",1), 
                       to = sapply(unlist(lapply(names(zw),strsplit,":"),recursive = FALSE),"[[",2),
                       weights = as.numeric(zw),stringsAsFactors = FALSE)
      ge <- graph.edgelist(as.matrix(zp[,1:2]),directed = TRUE)
      E(ge)$weights <- zp[,3]
      return(ge)
  }

  rescale <- function(x,newspan = c(-1,1)) {
      oldrange <- max(x) - min(x); newrange <- max(newspan) - min(newspan)
      scale <- oldrange/newrange; shift <- min(x) - (scale*-1)
      return((x - shift)/scale)
  }  

  # compute the similarity between scopes, and find hubs in each cover set compounds mst
  scd <- NULL
  cgraphlist <- NULL
  for (i in 1:length(mincoverset)) {
    iidx <- mincoverset[i]
    j <- i + 1
    while(j>i & j <= length(mincoverset)) {
      jidx <- mincoverset[j]
      z <- data.frame(from = names(scopes)[iidx],
                      to = names(scopes)[jidx],
                      weights = length(intersect(scopes[[iidx]],scopes[[jidx]])) )
      if (is.null(scd)) scd <- z
      else scd <- rbind(scd,z)
      j <- j + 1
    }
    ge <- getCompoundScopeNetwork(paste(loc,names(scopes)[iidx],".RData",sep = ""))
    w <- E(ge)$weights; w <- -log(rank(w)/length(w))
    mst <- minimum.spanning.tree(ge,weights = w)
    # hubs are top 10 in degree
    if (vcount(mst) > 10) {
          ng <- sort(degree(mst),decreasing = TRUE,index.return = TRUE)   
          cgraphlist[[i]] <- V(mst)$name[ng$ix[1:10]]
      } else {
        cgraphlist[[i]] <- V(mst)$name[-which(V(mst)$name == "s")] 
      }
      names(cgraphlist)[i] <- names(scopes)[iidx]
    }

  #set graph up
  mst <- graph.edgelist(as.matrix(scd[,1:2]),directed = FALSE)
  E(mst)$weights <- scd$weights
  V(mst)$ctype <- rep("compound",vcount(mst))

  # add intermediate hub reactions
  newnodes <- unique(as.character(unlist(cgraphlist)))
  inodes <- as.character(unlist(cgraphlist))
  mst <- add.vertices(mst,nv = length(newnodes),attr = list(name = newnodes,ctype = rep("hub",length(newnodes)) ))
  nreps <- as.character(unlist(mapply(rep,names(cgraphlist),sapply(cgraphlist,length))))
  ed <- cbind(nreps,inodes)
  ed <- as.character(t(unlist(ed)))
  ed <- match(ed,V(mst)$name)-1
  mst <- add.edges(mst,edges = ed,attr = list(weights = rep(10,length(inodes))))

  #node properties
  vcol <- rep("green",vcount(mst))
  vcol[V(mst)$ctype == "reaction"] <- "red"
  vcol[V(mst)$ctype == "hub"] <- "orange"
  vname <- rep("",vcount(mst))
  vname[V(mst)$ctype %in% "compound"] <- V(mst)$name[V(mst)$ctype %in% "compound"]
  vname[V(mst)$ctype %in% "hub"] <- V(mst)$name[V(mst)$ctype %in% "hub"]
  vshape <- "circle"
  vlab.col <- ifelse(V(mst)$ctype %in% "compound","blue","black")
  vsize <- rep(2,vcount(mst))
  ssize <- sapply(scopes[mincoverset],length)
  vsize[V(mst)$ctype %in% "compound"] <- ssize/max(ssize)*5 + 3
  hsize <- degree(mst)[V(mst)$ctype %in% "hub"]
  vsize[V(mst)$ctype %in% "hub"] <- hsize/max(hsize)*3 + 2
  
  #edge properties
  ecol <- rep("green",ecount(mst))
  ed <- suppressWarnings(get.edgelist(mst))
  ecol[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"hub"] <- "orange"
  ecol[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"reaction"] <- "red"
  ewidth <- rep(2,ecount(mst))
  ewidth[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"compound"] <- 6
  ewidth[V(mst)$ctype[match(ed[,2],V(mst)$name)]%in%"hub"] <- 3

  E(mst)$col <- ecol; V(mst)$col = vcol; V(mst)$shape <- vshape;V(mst)$vlab.col <- vlab.col
  V(mst)$size <- vsize; V(mst)$vname <- vname
  mst <- delete.edges(mst,E(mst)[ecol == "green" | ecol == "red"])
  mst <- delete.vertices(mst,V(mst)[vcol == "red"])
  mst <- delete.vertices(mst,V(mst)[degree(mst) == 0])
  
  zc <- clusters(mst)
  blist <- list(); idx <- 1
  cidx <- order(zc$csize,decreasing = TRUE)
  for (i in cidx) {
    gec <- subgraph(mst,V(mst)[zc$membership == (i-1)])
    el <- data.frame(get.edgelist(gec),stringsAsFactors = FALSE)
    names(el) <- c("from","to")
    blist[[idx]] <- bipartiteLayout(el,iter = 20)
    idx <- idx + 1
  }
  
  #layout the graphs
  zz <- zc$csize[cidx]/vcount(mst)
  zz[which.max(zz)] <- zz[which.max(zz)] - sum(zz < 0.05)*0.05
  zz[zz < 0.05] <- 0.05; zz <- zz/sum(zz)
  z <- layout(matrix(1:(length(blist)*3),length(blist),3,byrow = TRUE),widths = c(0.4,0.2,0.4),heights = zz)

  node.size <- degree(mst)/max(degree(mst)) * 5

  par(oma = c(2,2,3,2))
  for (i in 1:length(blist)) {
    zg <- blist[[i]]
    vc.idx <- match(zg$n1.pos[,1],V(mst)$name)
    vh.idx <- match(zg$n2.pos[,1],V(mst)$name)
    mn <- max(nrow(zg$n1.pos),nrow(zg$n2.pos))    
    yl <- c(-1 - 1/mn,1+1/mn)

    #Compound Labels
    par(mar = c(1,0,0,0))
    plot(NA,ylim = yl,xlim =c(0,1),axes = FALSE,xlab = "",ylab = "",yaxs = "i")
    vvn <- strsplit(sbml.info$species.info[match(V(mst)$name[vc.idx],sbml.info$species.info$id),]$name,"_")
    vidx <- strsplit(sbml.info$species.info[match(V(mst)$name[vc.idx],sbml.info$species.info$id),]$id,"_")
    vn <- c()
    for (k in 1:length(vvn)) {
      vcpd <- vidx[[k]][vidx[[k]] != ""]
      vcpd <- vcpd[length(vcpd)]
      z <- paste(vvn[[k]][vvn[[k]] != ""][1],"_",vcpd,sep = "")
      cid <- sbml.info$species.info[match(V(mst)$name[vc.idx[k]],sbml.info$species.info$id),]$id
      if (nchar(z) > 40) z <- paste(substr(z,0,40-nchar(cid)-5),"...(",cid,")",sep ="")
      vn <- c(vn,z)
    }    
    text(x = rep(1,nrow(zg$n1.pos)), y = zg$n1.pos$position,labels = vn,pos= 2,cex = 1.5)
    if (i == 1) mtext(side = 3,line = 0.5,"Cover Set Metabolites",cex =1.5,font =2)

    # draw the graph    
    plot(NA,xlim = c(-0.2,1.2),ylim = yl,axes = FALSE,xlab = "",ylab = "",yaxs = "i")
    
    # Position compounds
    nscl <- V(mst)$size[vc.idx]
    points(rep(0,nrow(zg$n1.pos)),y = zg$n1.pos$position,pch = 22,cex = nscl,
           col = "green",bg = "green")
    
    # Position reactions
    nscl <- V(mst)$size[vh.idx]
    points(rep(1,nrow(zg$n2.pos)),y = zg$n2.pos$position,pch = 21,
           cex = nscl,col = "orange",bg = "orange")
    
    # Edges
    ew <- rep(1,nrow(zg$edge.map))
    ieds <- unique(zg$edge.map[,2])
    if (length(ieds) > 2) {
      for (ee in ieds) ew[zg$edge.map[,2]==ee] <- sum(zg$edge.map[,2] == ee)
      ew <- ew/max(ew)*2
      ew[ew < 0.1] <- 0.1
    }
    segments(x0 = rep(0,nrow(zg$edge.map)),y0 = zg$edge.map[,1],
             x1 = rep(1,nrow(zg$edge.map)),y1 = zg$edge.map[,2],col = 2,lwd = ew)

    # Reaction Labels
    plot(NA,ylim = yl,xlim =c(0,1),axes = FALSE,xlab = "",ylab = "",yaxs = "i")
    vvn <- strsplit(sapply(sbml.info$reaction.info[V(mst)$name[vh.idx]],"[[","name"),"_")
    vn <- c()
    for (k in 1:length(vvn)) {
      z <- vvn[[k]][vvn[[k]] != ""][1]
      cid <- names(vvn)[k]
      if (nchar(z) > 40) z <- paste(substr(z,0,40-nchar(cid)-5),"...(",cid,")",sep ="")
      vn <- c(vn,z)
    }    
    text(x = rep(0,nrow(zg$n2.pos)), y = zg$n2.pos$position,labels = vn,pos= 4,cex = 1.5)
    if (i == 1) mtext(side = 3,line = 0.5,"Hub Reactions",cex =1.5,font =2)
  }

  return(mst)
}

