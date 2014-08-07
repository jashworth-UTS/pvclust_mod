hc2axes <- function(x)
{
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height

  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)

    for(i in 1:(n-1)) {
        ai <- A[i,1]

        if(ai < 0)
          x.tmp[1] <- zz[-ai]
        else
          x.tmp[1] <- x.axis[ai]

        ai <- A[i,2]

        if(ai < 0)
          x.tmp[2] <- zz[-ai]
        else
          x.tmp[2] <- x.axis[ai]

        x.axis[i] <- mean(x.tmp)
      }

  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

# for improved speed
hc2split = function(hc,return_members=T){
	.Call(hc2_split, hc$merge, return_members)
}

# deprecated/replaced by much faster c++ function above
hc2split_R <- function(x)
  {
	# VERY SLOW on big data; nearly instantaneous when rewritten in c++
	cat('hc2split...')
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()

    for(i in 1:(n-1)){
        ai <- A[i,1]

        if(ai < 0)
          B[[i]] <- -ai
        else
          B[[i]] <- B[[ai]]

        ai <- A[i,2]

        if(ai < 0)
          B[[i]] <- sort(c(B[[i]],-ai))
        else
          B[[i]] <- sort(c(B[[i]],B[[ai]]))
      }

    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)

    for(i in 1:(n-1)){
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
          CC[i,bi[j]] <- 1
      }

    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)

	  cat('done\n')
    return(split)
  }

pvclust.node <- function(x, r,...)
  {
#    require(pvclust)
    mboot.node <- lapply(r, boot.hclust, nboot=x, ...)
    return(mboot.node)
  }

boot.hclust_R <- function(r, data, object.hclust, method.dist, use.cor,
                        method.hclust, nboot, store, weight=F)
{
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n

	# following lines are TERRIBLY SLOW on big data. Very inefficient
	# factor(pattern) is very slow for sorting of large 'pattern' strings in R
	cat('init edges.cnt...')
  pattern   <- hc2split(object.hclust)$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
	cat('done\n')

  st <- list()

  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0

  for(i in 1:nboot){
		cat('.')
    if(weight && r>10) {  ## <- this part should be improved
      w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
      suppressWarnings(distance <- distw.pvclust(data,w1,method=method.dist,use.cor=use.cor))
    } else {
      smpl <- sample(1:n, size, replace=TRUE)
      suppressWarnings(distance  <- dist.pvclust(data[smpl,],method=method.dist,use.cor=use.cor))
    }
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust)
      pattern.i <- hc2split(x.hclust)$pattern # split
			# the following is terribly slow for large trees (10,000+)
			cat('counts...')
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
			cat('done\n')
    } else {
      x.hclust <- NULL
	  na.flag <- 1
    }

    if(store)
      st[[i]] <- x.hclust
  }
  cat("Done.\n")
  # bootstrap done

  if(na.flag == 1)
	warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)

  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"

  return(boot)
}

# modifed for major efficiency(speed) improvements in c++
boot.hclust <- function(r, data, object.hclust, method.dist, use.cor,
                        method.hclust, nboot, store, weight=F)
{
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n

  hcsplit = hc2split(object.hclust)

	# just keep track of unsorted tree correspondence (new cpp function below) and sort once later
	edges.cnt = rep(0,length(hcsplit$member))
	names(edges.cnt) = hcsplit$pattern

  st <- list()

  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0

  for(i in 1:nboot){
		cat('.')
    if(weight && r>10) {  ## <- this part should be improved
      w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
      distance <- distw.pvclust(data,w1,method=method.dist,use.cor=use.cor)
    } else {
      smpl <- sample(1:n, size, replace=TRUE)
      distance  <- dist.pvclust(data[smpl,],method=method.dist,use.cor=use.cor)
    }
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust <- hclust(distance,method=method.hclust)
			# keep count by comparing 'member' structures directly (c++ call)
			newcounts = .Call(count_edge_matches, object.hclust$merge, x.hclust$merge)
			oldec = edges.cnt
			edges.cnt = edges.cnt + newcounts
    } else {
      x.hclust <- NULL
			na.flag <- 1
    }

    if(store) st[[i]] <- x.hclust
  }
  cat("Done.\n")
  # bootstrap done

  if(na.flag == 1)
	warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)

	# now sort edges.cnt by 'pattern' string
	edges.cnt = edges.cnt[ order(names(edges.cnt)) ]

  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"

  return(boot)
}

pvclust.merge <- function(data, object.hclust, mboot){

  pattern <- hc2split(object.hclust)$pattern

  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")

  rl <- length(mboot)
  ne <- length(pattern)

  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")

  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt)
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }

  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"

  p    <- lapply(ms.fitted,"[[","p")
  se   <- lapply(ms.fitted,"[[","se")
  coef <- lapply(ms.fitted,"[[","coef")

  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))

  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)

  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne

  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store)

  class(result) <- "pvclust"
  return(result)
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
#	cat('Distances...\n')
  if(!is.na(pmatch(method,"correlation"))){
		#ja here substituting with a faster c++ method, which also more literally obeys the notion of 'pairwise.complete.obs' for Pearson correlation than the R source does
		require(fastcluster)
		res = fastcluster_correlation_distance(x)
#    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  }
  else
    dist(t(x),method)
}


corw <- function(x,w,
                 use=c("all.obs","complete.obs","pairwise.complete.obs")
                 ) {
  if(is.data.frame(x)) x <- as.matrix(x)
  x <- x[w>0,,drop=F]
  w <- w[w>0]

  n <- nrow(x) # sample size
  m <- ncol(x) # number of variables
  if(missing(w)) w <- rep(1,n)
  r <- matrix(0,m,m,dimnames=list(colnames(x),colnames(x)))
  diag(r) <- 1
  use <- match.arg(use)

  pairu <- F
  if(use=="all.obs") {
    u <- rep(T,n)
  } else if(use=="complete.obs") {
    u <- apply(x,1,function(y) !any(is.na(y)))
  } else if(use=="pairwise.complete.obs") {
    pairu <- T
    ux <- is.finite(x)
  } else stop("unknown use")

  for(i in 1+seq(length=m-1)) {
    for(j in seq(length=i-1)) {
      if(pairu) u <- ux[,i] & ux[,j]
      wu <- w[u]; xi <- x[u,i]; xj <- x[u,j]
      ws <- sum(wu)
      if(ws > 1e-8) {
        xi <- xi - sum(wu*xi)/ws
        xj <- xj - sum(wu*xj)/ws
        vxi <- sum(wu*xi*xi)/ws
        vxj <- sum(wu*xj*xj)/ws
        if(min(vxi,vxj) > 1e-8)  {
          vxij <- sum(wu*xi*xj)/ws
          rij <- vxij/sqrt(vxi*vxj)
        } else {
          rij <- 0
        }
      } else {
        rij <- 0
      }
      r[i,j] <- r[j,i] <- rij
    }
  }
  r
}

### calculate distance by weight
distw.pvclust <- function(x,w,method="correlation", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - corw(x,w, use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(corw(x,w, use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  stop("wrong method")
}
