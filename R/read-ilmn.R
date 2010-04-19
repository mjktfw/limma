.readData <- function(fname, probeid, annotation, expr, other.columns, sep, quote, verbose, ...)
{
skiplines <- 0
repeat{
  rl <- scan(fname, what="", sep=sep, quote=quote, nlines=1, quiet=TRUE, skip=skiplines, ...)
  if(length(grep(tolower(expr), tolower(rl)))) break
  skiplines <- skiplines + 1 
}
header <- unlist(strsplit(rl, sep))

elist <- new("EListRaw")
reqcol <- header[grep(tolower(paste(c(probeid, annotation, expr, other.columns), collapse="|")), tolower(header))]

x <- read.columns(file=fname, required.col=reqcol, skip=skiplines, sep=sep, quote=quote, stringsAsFactors=FALSE,  ...)
nprobes <- nrow(x)

pids <- x[, grep(tolower(probeid), tolower(colnames(x)))]
snames <- colnames(x)[grep(tolower(expr), tolower(colnames(x)))]
snames <- unlist(strsplit(snames, paste("[.]*", expr, "-*", sep="")))
snames <- snames[snames != ""]
nsamples <- length(snames)

elist$E <- data.matrix(x[, grep(tolower(expr), tolower(colnames(x)))])
colnames(elist$E) <- snames
rownames(elist$E) <- pids

elist$genes <- x[, c(grep(tolower(probeid), tolower(colnames(x))), grep(tolower(paste(annotation,collapse="|")), tolower(colnames(x))))]
elist$targets <- data.frame(SampleNames=snames, stringsAsFactors=FALSE)

if(!is.null(other.columns)){
  elist$other <- vector("list", length(other.columns))
  for(i in 1:length(other.columns)){
    if(length(co <- grep(tolower(other.columns[i]), tolower(colnames(x)))))
      elist$other[[i]] <- as.matrix(x[, co])
    else
      elist$other[[i]] <- matrix(NA, nprobes, nsamples)
    colnames(elist$other[[i]]) <- snames
    rownames(elist$other[[i]]) <- pids
  }
  names(elist$other) <- other.columns
}
elist
}

read.ilmn <- function(files=NULL, ctrlfiles=NULL, path=NULL, ctrlpath=NULL,
probeid="Probe",  annotation=c("TargetID", "SYMBOL"), expr="AVG_Signal", other.columns=NULL, 
sep="\t", quote="\"", verbose=TRUE, ...)
{
if(is.null(ctrlfiles))
  warning("Names of control profile files are not provided.")

if(!is.null(files)){
  f <- unique(files)
  if(!is.null(path)) f <- file.path(path, f)
  n <- length(f)
  for(i in 1:n){
    if(verbose) cat("Reading file", f[i], "... ...\n") 
    elist1 <- .readData(f[i], probeid,  annotation, expr, other.columns, sep, quote, verbose, ...)
    if(i==1)
      elist <- elist1
    else
      elist <- cbind(elist, elist1)
  }
  if(!is.null(ctrlfiles))
    elist$genes$Status <- "regular"
}
  
if(!is.null(ctrlfiles)){
  cf <- unique(ctrlfiles)
  if(!is.null(ctrlpath)) cf <- file.path(ctrlpath, cf)
  n <- length(cf)
  for(i in 1:n){
    if(verbose) cat("Reading file", cf[i], "... ...\n")
    elist.ctrl1 <- .readData(cf[i], probeid,  annotation, expr, other.columns, sep, quote, verbose, ...)
    if(i==1)
      elist.ctrl <- elist.ctrl1
    else
      elist.ctrl <- cbind(elist.ctrl, elist.ctrl1)
  }
  elist.ctrl$genes$Status <- elist.ctrl$genes[,ncol(elist.ctrl$genes)]
}

if(!is.null(files))
  if(!is.null(ctrlfiles)){
    colnames(elist.ctrl$genes) <- colnames(elist$genes)
    return(rbind(elist, elist.ctrl))
  }
  else
    return(elist)
else
  if(!is.null(ctrlfiles))
    return(elist.ctrl)
}

read.ilmn.targets <- function(targets, ...)
{
if(is.null(targets$files) && is.null(targets$ctrlfiles))
  stop("Can not find column \"files\" or \"ctrlfiles\" in targets\n") 
x <- read.ilmn(targets$files, targets$ctrlfiles, ...) 
x$targets <- targets
x
}

