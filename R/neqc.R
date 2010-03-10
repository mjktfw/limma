normexp.fit.control <-
function(x, status=NULL, negctrl="negative"){
if(is(x, "EListRaw")){
  if(!is.null(status))
    s <- status
  else
    if(is.null(s <- x$genes$Status))
      stop("Probe status can not be found!")
  xr <- x$E[tolower(s)=="regular", ]
  xn <- x$E[tolower(s)==tolower(negctrl), ]
}
else
  if(is.matrix(x)){
    if(is.null(status))
      stop("Please provide probe status!") 
    xr <- x[tolower(status)=="regular", ]
    xn <- x[tolower(status)==tolower(negctrl), ]
  }
  else
    stop("Data type not supported!\n")

mu <- sigma <- alpha <-  rep(NA, ncol(x))
for(i in 1:ncol(x)){
  mu[i] <- mean(xn[, i], na.rm=TRUE)
  sigma[i] <- sd(xn[, i], na.rm=TRUE)
  alpha[i] <- max(mean(xr[, i], na.rm=TRUE) - mu[i], 10)
}
normexp.par <- matrix(cbind(mu, log(sigma), log(alpha)), ncol=3)
colnames(normexp.par) <- c("mu", "logsigma", "logalpha")
return(normexp.par)
}

neqc <-
function(x, status=NULL, negctrl="negative", regular="regular", offset=16, ...){
normexp.par <- normexp.fit.control(x, status, negctrl)
if(is(x, "EListRaw")){
  for(i in 1:ncol(x))
    x$E[, i] <- normexp.signal(normexp.par[i, ], x$E[, i])
  x$E <- x$E + offset
  y <- normalizeBetweenArrays(x, method="quantile", ...)
  if(is.null(status))
    status <- x$genes$Status
  y <- y[status == regular, ]
  y$genes$Status <- NULL
}
if(is.matrix(x)){
  for(i in 1:ncol(x))
    x[, i] <- normexp.signal(normexp.par[i, ], x[, i])
  x <- x + offset
  y <- log2(normalizeBetweenArrays(x, method="quantile", ...))
  y <- y[status == regular, ]
}
y
}
