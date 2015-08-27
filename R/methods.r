#' @title initialize
#' @docType methods
#' @rdname initialize
setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )
#' @title filterdata
#' @export
#' @docType methods
#' @rdname filterdata
setGeneric("filterdata", function(object, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000) standardGeneric("filterdata"))

#' @describeIn filterdata
setMethod("filterdata",
          signature = "SCseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
            if ( downsample ){
              set.seed(rseed)
              object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
            }else{
              x <- object@ndata
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata
            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
            return(object)
          }
          )

#' @title clustexp
#' @export
#' @docType methods
#' @rdname clustexp
setGeneric("clustexp", function(object,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000) standardGeneric("clustexp"))

#' @describeIn clustexp
setMethod("clustexp",
          signature = "SCseq",
          definition = function(object,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed) {
            if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
            if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
            if ( ! ( metric %in% c( "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski") ) ) stop("metric has to be one of the following: spearman, pearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski")
            if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
            if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
            if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
            if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")          
            if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
            if ( !do.gap & cln == 0 ) stop("cln has to be a positive integer or do.gap has to be TRUE")
            object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed)
            y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed)
            object@kmeans    <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr)
            object@distances <- as.matrix( y$di )
            set.seed(111111)
            object@fcol <- sample(rainbow(max(y$clb$result$partition)))
            return(object)
          }
          )

#' @export
#' @docType methods
#' @rdname findoutliers
setGeneric("findoutliers", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75) standardGeneric("findoutliers"))

#' @describeIn findoutliers
setMethod("findoutliers",
          signature = "SCseq",
          definition = function(object,outminc,outlg,probthr,thr,outdistquant) {
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before findoutliers")
            if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
            if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
            if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
            if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
            if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
                      
            object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
            ### calibrate background model
            m <- log2(apply(object@fdata,1,mean))
            v <- log2(apply(object@fdata,1,var))
            f <- m > -Inf & v > -Inf
            m <- m[f]
            v <- v[f]
            mm <- -8
            repeat{
              fit <- lm(v ~ m + I(m^2)) 
              if( coef(fit)[3] >= 0 | mm >= 3){
                break
              }
              mm <- mm + .5
              f <- m > mm
              m <- m[f]
              v <- v[f]
            }
            object@background <- list()
            object@background$vfit <- fit
            object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
            object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)

            ### identify outliers
            out   <- c()
            stest <- rep(0,length(thr))
            cprobs <- c()
            for ( n in 1:max(object@kmeans$kpart) ){     
              if ( sum(object@kmeans$kpart == n) == 1 ){ cprobs <- append(cprobs,.5); names(cprobs)[length(cprobs)] <- names(object@kmeans$kpart)[object@kmeans$kpart == n]; next }
              x <- object@fdata[,object@kmeans$kpart == n]
              x <- x[apply(x,1,max) > outminc,]
              z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
              cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
              f <- cp < probthr
              cprobs <- append(cprobs,cp)
              if ( sum(f) > 0 ) out <- append(out,names(x)[f])
              for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )  
            }
            object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs)

            ### cluster outliers
            clp2p.cl <- c()
            cols <- names(object@fdata)
            di <- as.data.frame(object@distances)
            for ( i in 1:max(object@kmeans$kpart) ) {
              tcol <- cols[object@kmeans$kpart == i]
              if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
            }
            clp2p.cl <- clp2p.cl[clp2p.cl>0]
  
            cpart <- object@kmeans$kpart
            cadd  <- list()
            if ( length(out) > 0 ){
              if (length(out) == 1){
                cadd <- list(out)
              }else{
                n <- out
                m <- as.data.frame(di[out,out])
                
                for ( i in 1:length(out) ){
                  if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]          
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1  
                    if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
          
                  }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                  }
                }
              }
    
              for ( i in 1:length(cadd) ){
                cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
              }
            }

            ### determine final clusters
            for ( i in 1:max(cpart) ){
              d <- object@fdata[,cols[cpart == i]]
              if ( sum(cpart == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
              if ( i == 1 ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
              if ( i == 1 ) tmp <- data.frame(apply(object@fdata[,cols],2,dist.gen.pairs,y=cent,method=object@clusterpar$metric)) else tmp <- cbind(tmp,apply(object@fdata[,cols],2,dist.gen.pairs,y=cent,method=object@clusterpar$metric))
            }    
            cpart <- apply(tmp,1,function(x) order(x,decreasing=FALSE)[1])
  
            for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }

            object@cpart <- cpart

            set.seed(111111)
            object@fcol <- sample(rainbow(max(cpart)))
            return(object)
          }
        )

#' @export
#' @docType methods
#' @rdname plotgap
setGeneric("plotgap", function(object) standardGeneric("plotgap"))

#' @describeIn plotgap
setMethod("plotgap",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@kmeans$gap)
          }
          )

#' @export
#' @docType methods
#' @rdname plotsilhouette
setGeneric("plotsilhouette", function(object) standardGeneric("plotsilhouette"))

#' @describeIn plotsilhouette
setMethod("plotsilhouette",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@kmeans$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            kpart <- object@kmeans$kpart
            distances  <- dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si)
          }
          )

#' @export
#' @docType methods
#' @rdname plotjaccard
setGeneric("plotjaccard", function(object) standardGeneric("plotjaccard"))

#' @describeIn plotjaccard
setMethod("plotjaccard",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotjaccard")
            if ( length(unique(object@kmeans$kpart)) < 2 ) stop("only a single cluster: no Jaccard's similarity plot")
            barplot(object@kmeans$jaccard,names.arg=1:length(object@kmeans$jaccard),ylab="Jaccard's similarity")
          }
          )

#' @export
#' @docType methods
#' @rdname plotoutlierprobs
setGeneric("plotoutlierprobs", function(object) standardGeneric("plotoutlierprobs"))

#' @describeIn plotoutlierprobs
setMethod("plotoutlierprobs",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotoutlierprobs")
            p <- object@kmeans$kpart[ order(object@kmeans$kpart,decreasing=FALSE)]
            x <- object@out$cprobs[names(p)]
            fcol <- object@fcol
            for ( i in 1:max(p) ){
              y <- -log10(x + 2.2e-16)
              y[p != i] <- 0
              if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
  }
            abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
            d <- b[2,1] - b[1,1]
            y <- 0
            for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,lab=1:max(p))
            box()
          }
          )

#' @export
#' @docType methods
#' @rdname plotbackground
setGeneric("plotbackground", function(object) standardGeneric("plotbackground"))

#' @describeIn plotbackground
setMethod("plotbackground",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotbackground")
            m <- apply(object@fdata,1,mean)
            v <- apply(object@fdata,1,var)
            fit <- locfit(v~lp(m,nn=.7),family="gamma",maxk=500)
            plot(log2(m),log2(v),pch=20,xlab="log2mean",ylab="log2var",col="grey")
            lines(log2(m[order(m)]),log2(object@background$lvar(m[order(m)],object)),col="red",lwd=2)
            lines(log2(m[order(m)]),log2(fitted(fit)[order(m)]),col="orange",lwd=2,lty=2)
            legend("topleft",legend=substitute(paste("y = ",a,"*x^2 + ",b,"*x + ",d,sep=""),list(a=round(coef(object@background$vfit)[3],2),b=round(coef(object@background$vfit)[2],2),d=round(coef(object@background$vfit)[3],2))),bty="n")
          }
          )

#' @export
#' @docType methods
#' @rdname plotsensitivity
setGeneric("plotsensitivity", function(object) standardGeneric("plotsensitivity"))

#' @describeIn plotsensitivity
setMethod("plotsensitivity",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotsensitivity")
            plot(log10(object@out$thr), object@out$stest, type="l",xlab="log10 Probability cutoff", ylab="Number of outliers")
            abline(v=log10(object@outlierpar$probthr),col="red",lty=2)
          }
          )

#' @title clustdiffgenes
#' @name clustdiffgenes
#' @export
#' @docType methods
#' @rdname clustdiffgenes
setGeneric("clustdiffgenes", function(object,pvalue=.01) standardGeneric("clustdiffgenes"))

#' @describeIn clustdiffgenes
setMethod("clustdiffgenes",
          signature = "SCseq",
          definition = function(object,pvalue){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
            if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@cpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <- apply(x,1,mean)
              n <- if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
            }
            return(cdiff)
          }
          )

#' @export
#' @docType methods
#' @rdname diffgenes
setGeneric("diffgenes", function(object,cl1,cl2,mincount=5) standardGeneric("diffgenes"))

#' @describeIn diffgenes
setMethod("diffgenes",
          signature = "SCseq",
          definition = function(object,cl1,cl2,mincount){
            part <- object@cpart
            cl1 <- c(cl1)
            cl2 <- c(cl2)
            if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
            if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
            if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            f <- apply(object@ndata[,part %in% c(cl1,cl2)],1,max) > mincount
            x <- object@ndata[f,part %in% cl1]
            y <- object@ndata[f,part %in% cl2]
            if ( sum(part %in% cl1) == 1 ) m1 <- x else m1 <- apply(x,1,mean)
            if ( sum(part %in% cl2) == 1 ) m2 <- y else m2 <- apply(y,1,mean)
            if ( sum(part %in% cl1) == 1 ) s1 <- sqrt(x) else s1 <- sqrt(apply(x,1,var))
            if ( sum(part %in% cl2) == 1 ) s2 <- sqrt(y) else s2 <- sqrt(apply(y,1,var))
            
            d <- ( m1 - m2 )/ apply( cbind( s1, s2 ),1,mean )
            names(d) <- rownames(object@ndata)[f]
            if ( sum(part %in% cl1) == 1 ){
              names(x) <- names(d)
              x <- x[order(d,decreasing=TRUE)]
            }else{
              x <- x[order(d,decreasing=TRUE),]
            }
            if ( sum(part %in% cl2) == 1 ){
              names(y) <- names(d)
              y <- y[order(d,decreasing=TRUE)]
            }else{
              y <- y[order(d,decreasing=TRUE),]
            }
            return(list(z=d[order(d,decreasing=TRUE)],cl1=x,cl2=y,cl1n=cl1,cl2n=cl2))
          }
          )

#' @export
#' @docType methods
#' @rdname clustheatmap
setGeneric("clustheatmap", function(object,final=FALSE,hmethod="single") standardGeneric("clustheatmap"))

#' @describeIn clustheatmap
setMethod("clustheatmap",
          signature = "SCseq",
          definition = function(object,final,hmethod){
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before clustheatmap")
            if ( !final & length(object@kmeans$kpart) == 0 ) stop("run clustexp before clustheatmap")
            x <- object@fdata  
            part <- if ( final ) object@cpart else object@kmeans$kpart
            na <- c()
            j <- 0
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              j <- j + 1
              na <- append(na,i)
              d <- x[,part == i]
              if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
              if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
            }
            names(tmp) <- paste("cl",na,sep=".")
            if ( max(part) > 1 )  cclmo <- hclust(dist.gen(as.matrix(dist.gen(t(tmp),method=object@clusterpar$metric))),method=hmethod)$order else cclmo <- 1
            q <- part
            for ( i in 1:max(part) ){
              q[part == na[cclmo[i]]] <- i
            }
            part <- q
            di <- as.data.frame( as.matrix( dist.gen(t(object@distances)) ) )
            pto <- part[order(part,decreasing=FALSE)]
            ptn <- c()
            for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
            col <- object@fcol
            mi  <- min(di,na.rm=TRUE)
            ma  <- max(di,na.rm=TRUE)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            if ( mi == ma ){
              ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
            }
            par(mar = c(3,5,2.5,2))
            image(as.matrix(di[ptn,ptn]),col=ColorRamp,axes=FALSE)
            abline(0,1)
            box()
            
            tmp <- c()
            for ( u in 1:max(part) ){
              ol <- (0:(length(part) - 1)/(length(part) - 1))[ptn %in% names(x)[part == u]]
              points(rep(0,length(ol)),ol,col=col[cclmo[u]],pch=15,cex=.75)
              points(ol,rep(0,length(ol)),col=col[cclmo[u]],pch=15,cex=.75)
              tmp <- append(tmp,mean(ol))
            }
            axis(1,at=tmp,lab=cclmo)
            axis(2,at=tmp,lab=cclmo)
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
            return(cclmo)
          }
          )

#' @export
#' @docType methods
#' @rdname comptsne
setGeneric("comptsne", function(object,rseed=15555) standardGeneric("comptsne"))

#' @describeIn comptsne
setMethod("comptsne",
          signature = "SCseq",
          definition = function(object,rseed){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before comptsne")
            set.seed(rseed)
            di <- dist.gen(as.matrix(object@distances))
            ts <- tsne(di,k=2)
            object@tsne <- as.data.frame(ts)
            return(object)
          }
          )

#' @export
#' @docType methods
#' @rdname plottsne
setGeneric("plottsne", function(object,final=TRUE) standardGeneric("plottsne"))

#' @describeIn plottsne
setMethod("plottsne",
          signature = "SCseq",
          definition = function(object,final){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
            if ( !final & length(object@kmeans$kpart) == 0 ) stop("run clustexp before plottsne")
            part <- if ( final ) object@cpart else object@kmeans$kpart
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=object@fcol[i],cex=.75,font=4)
            }
          }
          )

#' @export
#' @docType methods
#' @rdname plotlabelstsne
setGeneric("plotlabelstsne", function(object,labels=NULL) standardGeneric("plotlabelstsne"))

#' @describeIn plotlabelstsne
setMethod("plotlabelstsne",
          signature = "SCseq",
          definition = function(object,labels){
            if ( is.null(labels ) ) labels <- names(object@ndata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotlabelstsne")
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            text(object@tsne[,1],object@tsne[,2],labels,cex=.5)
          }
          )

#' @export
#' @docType methods
#' @rdname plotsymbolstsne
setGeneric("plotsymbolstsne", function(object,types=NULL) standardGeneric("plotsymbolstsne"))

#' @describeIn plotsymbolstsne
setMethod("plotsymbolstsne",
          signature = "SCseq",
          definition = function(object,types){
            if ( is.null(types) ) types <- names(object@fdata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotsymbolstsne")
            if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
            coloc <- rainbow(length(unique(types)))
            syms <- c()
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,col="grey")
            for ( i in 1:length(unique(types)) ){
              f <- types == sort(unique(types))[i]
              syms <- append( syms, ( (i-1) %% 25 ) + 1 )
              points(object@tsne[f,1],object@tsne[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
            }
            legend("topleft", legend=sort(unique(types)), col=coloc, pch=syms)
          }
          )

#' @export
#' @docType methods
#' @rdname plotexptsne
setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE) standardGeneric("plotexptsne"))

#' @describeIn plotexptsne
setMethod("plotexptsne",
          signature = "SCseq",
          definition = function(object,g,n="",logsc=FALSE){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            if (logsc) {
              f <- l == 0
              l <- log(l)
              l[f] <- NA
            }
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey")
            for ( k in 1:length(v) ){
              points(object@tsne[k,1],object@tsne[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
            }
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
          }
          )
