#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="subtract") {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 9 Dec 2003.

	method <- match.arg(method, c("none","subtract", "half", "minimum", "edwards"))
	switch(method,
	subtract={
		RG$R <- RG$R-RG$Rb
		RG$G <- RG$G-RG$Gb
	},
	half={
		RG$R <- pmax(RG$R-RG$Rb, 0.5)
		RG$G <- pmax(RG$G-RG$Gb, 0.5)
	},
	minimum={
		RG$R <- as.matrix(RG$R - RG$Rb)
		RG$G <- as.matrix(RG$G - RG$Gb)
		for (slide in 1:ncol(RG$R)) {
			i <- RG$R[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$R[!i,slide],na.rm=TRUE)
				RG$R[i,slide] <- m/2
			}
			i <- RG$G[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$G[!i,slide],na.rm=TRUE)
				RG$G[i,slide] <- m/2
			}
		}
	},
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
		nspots <- NROW(RG$R)
		del <- function(d,f=0.1) quantile(d,sum(d<1e-16)*(1+f)/length(d),na.rm=TRUE)
		sub <- as.matrix(RG$R-RG$Rb)
		delta <- apply(sub, 2, del)
		RG$R <- ifelse(sub > rep(1,nspots)%o%delta, sub, delta*exp(1-(RG$Rb+delta)/RG$R))
		sub <- as.matrix(RG$G-RG$Gb)
		delta <- apply(sub, 2, del)
		RG$G <- ifelse(sub > rep(1,nspots)%o%delta, sub, delta*exp(1-(RG$Gb+delta)/RG$G))
	})
	RG$Rb <- NULL
	RG$Gb <- NULL
	new("RGList",unclass(RG))
}
