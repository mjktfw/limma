#  QUALITY MEASURES

wtarea <- function(ideal=c(160,170))
#	Quality weights based on spot area from SPOT output
#	Gordon Smyth
#	9 March 2003.  Last revised 11 Mar 2003.

function(spot) {
	e <- range(ideal)
	x <- c(-Inf,0,e,sum(e),Inf)
	y <- c(0,0,1,1,0,0)
	approx(x,y,xout=spot[,"area"],ties="ordered")$y
}

wtflags <- function(flagged.wt=0.1)
#	Quality weights based on Flags from GenePix output
#	Gordon Smyth
#	9 March 2003.  Last revised 18 Mar 2003.

function(gpr) {
	flagged <- (gpr[,"Flags"] < 0)
	flagged.wt*flagged + !flagged
}
