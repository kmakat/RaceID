#' @export
clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000)
{
  if ( clustnr < 2) stop("Choose clustnr > 1")
  di <- dist.gen(t(x),method=metric)
  if ( do.gap | cln > 0 ){
    gpr <- NULL
    if ( do.gap ){
      set.seed(rseed)
      gpr <- clusGap(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap) 
      if ( cln == 0 ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
    }    
    if ( cln <= 1 ) {
      clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
      names(clb$result$partition) <- names(x)
      return(list(x=x,clb=clb,gpr=gpr,di=di))
    }
    clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=KmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    return(list(x=x,clb=clb,gpr=gpr,di=di))
  }
}