contrastAsCoef <- function(design, contrast=NULL, first=TRUE)
#	Reform a design matrix so that one or more contrasts become simple coefficients
#	Gordon Smyth
#	31 August 2013
{
	design <- as.matrix(design)
	if(is.null(contrast)) return(design)
	contrast <- as.matrix(contrast)
	if(ncol(design) != nrow(contrast)) stop("Length of contrast doesn't match ncol(design)")
	qrc <- qr(contrast)
	ncontrasts <- qrc$rank
	if(ncontrasts==0) stop("contrast is all zero")
	coef <- 1:ncontrasts
	Dvec <- diag(qrc$qr)[coef]
	design <- t(qr.qty(qrc,t(design)))
	colnames(design) <- paste("Q",1:ncol(design),sep="")
	cn <- colnames(contrast)
	if(is.null(cn)) cn <- paste("C",qrc$pivot[coef],sep="")
	colnames(design)[coef] <- cn
	if(!first) {
		design <- cbind(design[,-coef,drop=FALSE],design[,coef,drop=FALSE])
		coef <- rev( ncol(design)-coef+1 )
	}
	design[,coef] <- t(t(design[,coef])/Dvec)
	list(design=design,coef=coef,qr=qrc)
}
