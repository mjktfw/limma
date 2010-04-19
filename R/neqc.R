normexp.fit.control <-
function(x, status=NULL, negctrl="negative", regular="regular", robust=FALSE)
#  Wei Shi
#  Edits by Gordon Smyth, 16 April 2010
#  Edits by Wei Shi, 19 April 2010
{
if(is(x, "EListRaw")){
  if(!is.null(status))
    s <- status
  else
    if(is.null(s <- x$genes$Status))
      stop("Probe status can not be found!")
  x$E <- as.matrix(x$E)
  xr <- x$E[tolower(s)==tolower(regular),,drop=FALSE]
  xn <- x$E[tolower(s)==tolower(negctrl),,drop=FALSE]
}
else {
  if(is.null(status)) stop("Please provide probe status!")
  x <- as.matrix(x)
  xr <- x[tolower(status)==tolower(regular),,drop=FALSE]
  xn <- x[tolower(status)==tolower(negctrl),drop=FALSE]
}

if(robust) {
	require(statmod)
	mu <- apply(xn,2,mean,trim=0.2,na.rm=TRUE)
	sigma <- apply(t(xn)-mu,1,mscale,na.rm=TRUE)
} else {
	mu <- colMeans(xn,na.rm=TRUE)
	sigma <- sqrt(rowSums((t(xn)-mu)^2,na.rm=TRUE)/(nrow(xn)-1))
}
alpha <- pmax(colMeans(xr,na.rm=TRUE)-mu,10)
cbind(mu=mu,logsigma=log(sigma),logalaph=log(alpha))
}

neqc <-
function(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE, ...)
#  Wei Shi
#  Edits by Gordon Smyth, 17 April 2010
{
normexp.par <- normexp.fit.control(x, status=status, negctrl=negctrl, regular=regular, robust=robust)
if(is(x, "EListRaw")) {
  for(i in 1:ncol(x))
    x$E[, i] <- normexp.signal(normexp.par[i, ], x$E[, i])
  x$E <- x$E + offset
  y <- normalizeBetweenArrays(x, method="quantile", ...)
  if(is.null(status))
    status <- x$genes$Status
  y <- y[tolower(status) == tolower(regular), ]
  y$genes$Status <- NULL
} else {
  x <- as.matrix(x)	
  for(i in 1:ncol(x))
    x[, i] <- normexp.signal(normexp.par[i, ], x[, i])
  x <- x + offset
  y <- log2(normalizeBetweenArrays(x, method="quantile", ...))
  y <- y[tolower(status) == tolower(regular), ]
}
y
}
