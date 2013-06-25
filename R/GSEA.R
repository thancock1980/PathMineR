plotEnrichmentScore <- function(eres) {
  gloc <- eres$gsea.res$geneset.locations

  layout(matrix(c(1,2,3),nrow =3,ncol =1),width = c(1,1,1),height = c(0.4,0.1,0.5))
  par(mar = c(0,4,4,4))
  zm <- apply(do.call("cbind",lapply(eres$perm.res,"[[","es")),1,max)
  mv <- ceiling(max(eres$gsea.res$es)*10)/10
  yl <- c( -1,1)
  plot(x = 1:length(gloc),y = eres$gsea.res$es,lwd = 2,col = 2,type = "l",ylab = "GSEA Enrichment Score",
       axes = FALSE,xaxs = "i",ylim = yl,yaxs = "i",main = "Enrichment Score Running Profile")
  axis(side = 2,at = seq(yl[1],yl[2],0.1),seq(yl[1],yl[2],0.1))
  lines(x = 1:length(gloc),y = zm,lwd = 1,col = 1)
  par(mar = c(4,4,0,4))
  image(x = 1:length(gloc),z = as.matrix(as.numeric(eres$gsea.res$geneset.locations)),xaxs = "i",
    ylab = "",xlab = "Selected Gene Pair Edges Sorted by Decreasing Positive Correlation",
    col = c("white","black"),axes = FALSE);box()
  axis(side = 1)

  par(mar = c(4,4,4,4))
  pnull <- sapply(eres$perm.res,"[[","enrichment.score")
  zh <- hist(pnull,main = "Enrichment Score Permuted Null Distribution",xlab = "Enrichment Score",xlim = c(0,1))
  abline(v=eres$gsea.res$enrichment.score,col = 2,lwd = 2)
  text(x = eres$gsea.res$enrichment.score,y = max(zh$counts),
      paste("P-value =",round(eres$pvalue,4)),pos = 2,col = 2)
}

enrichmentScore <- function(w,geneset,p=1,perms=1000) {
  hitmiss <- function(i,zw,isc,N,Nr,Nh,PERM = FALSE) {
    if (PERM) {  
        isc <- sample(isc)
        Nr <- sum(abs(zw[isc])^p) 
    }
   	idx <- which(isc)
	  zscore <- cumsum( abs(zw[which(isc)])^p )/Nr
	  hit <- unlist(mapply(rep,zscore[-length(zscore)],diff(idx)))
	  if (idx[1] != 1) hit <- c(rep(0,idx[1]-1),hit)
	  if (length(hit) != N) hit <- c(hit,rep(1,N-length(hit)))

	  idx <- which(!isc)
	  zscore <- cumsum(rep(1/(N-Nh),sum(!isc)))
	  miss <- unlist(mapply(rep,zscore[-length(zscore)],diff(idx)))
	  if (idx[1] != 1) miss <- c(rep(0,idx[1]-1),miss)
	  if (length(miss) != N) miss <- c(miss,rep(1,N-length(miss)))
	     
	  es <- hit - miss
    return(list(enrichment.score = max(es),es = es,hit = hit,miss = miss,geneset.locations = isc))
 }

 rew <- order(w,decreasing = TRUE)
 gs <- geneset[rew];wghts <- w[rew]
 zN <- length(rew);zNh <- sum(geneset);zNr <- sum(abs(wghts[gs])^p)
 gsea.res <- hitmiss(1,wghts,gs,zN,zNr,zNh,FALSE)
 perm.res <- lapply(1:perms,hitmiss,wghts,gs,zN,zNr,zNh,TRUE)
 pnull <- sapply(perm.res,"[[","enrichment.score")
 pnullecdf <- ecdf(pnull)
 pv <- pnullecdf(gsea.res$enrichment.score)

 return(list(gsea.res = gsea.res,perm.res = perm.res,pvalue = 1-pv))
}


