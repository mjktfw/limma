#  KOORPERBERG BACKGROUND ADUSTMENT FOR GENEPIX DATA

bayesianAdjustedRG <- function(slidenames, meanfg=TRUE, meanbg=FALSE)
	{
	RG <- list(R=NULL, G=NULL)
	for(i in 1:length(slidenames))
		{
		temp <- bayesianAdjustedFG(get(slidenames[i]), meanfg, meanbg)
		RG$R <- cbind(RG$R, temp$R)
		RG$G <- cbind(RG$G, temp$G)
		}
	cat(paste("Model-based background correction completed for", length(slidenames), "slides\n", sep=" "))
	RG
	}
	
bayesianAdjustedFG <- function(slide, meanfg=TRUE, meanbg=FALSE) {
	numspots <- dim(slide)[1]
	RG <- list(R=NULL, G=NULL)
        if(meanfg) {
          Rfg <- slide[,"F635.Mean"] 
          Gfg <- slide[,"F532.Mean"] 
        }
        if(!meanfg) {
          Rfg <- slide[,"F635.Median"] 
          Gfg <- slide[,"F532.Median"] 
        }
        if(meanbg) {
          Rbg <- slide[,"B635.Mean"] 
          Gbg <- slide[,"B532.Mean"] 
	  }  
	if(!meanbg) {
          Rbg <- slide[,"B635.Median"] 
          Gbg <- slide[,"B532.Median"] 
	  }
	Rsfg <- slide[,"F635.SD"]/sqrt(slide[,"F.Pixels"]) 
	Rsbg <- slide[,"B635.SD"]/sqrt(slide[,"B.Pixels"]) 
	Gsfg <- slide[,"F532.SD"]/sqrt(slide[,"F.Pixels"]) 
	Gsbg <- slide[,"B532.SD"]/sqrt(slide[,"B.Pixels"]) 
        for(i in 1:numspots) { # numspots
              RG$R <- c(RG$R, expectedBayesianAdjustedFG(fg=Rfg[i],
                                bg=Rbg[i], sfg=Rsfg[i], sbg=Rsbg[i]))
              RG$G <- c(RG$G, expectedBayesianAdjustedFG(fg=Gfg[i],
                                bg=Gbg[i], sfg=Gsfg[i], sbg=Gsbg[i]))
              }
        dim(RG$R) <- c(numspots,1)
        dim(RG$G) <- c(numspots,1)
        RG
      }
      
expectedBayesianAdjustedFG <- function(fg, bg, sfg, sbg) {
	integrate(numeratorBayesianAdjustedFG, ifelse((fg-bg-4*sqrt(sbg^2+sfg^2))<0, 0, fg-bg-4*sqrt(sbg^2+sfg^2)), 
		ifelse((fg-bg+4*sqrt(sfg^2+sbg^2))<0, 1000, fg-bg+4*sqrt(sfg^2+sbg^2)) , fg=fg, bg=bg, sfg=sfg, sbg=sbg, subdivisions=10000)$value/denominatorBayesianAdjustedFG(fg, bg, sfg, sbg)
	}
	
numeratorBayesianAdjustedFG <- function(ut, fg, bg, sfg, sbg)
	ut*exp(dnorm((fg-ut-bg)/sqrt(sfg^2+sbg^2), log=TRUE)+pnorm(((fg-ut)*sbg^2+bg*sfg^2)/(sbg*sfg*sqrt(sfg^2+sbg^2)), log.p=TRUE))
	
denominatorBayesianAdjustedFG <- function(fg, bg, sfg, sbg)	{
	normalConvolution <- function(v, fg, bg, sfg, sbg) exp(pnorm((fg-v)/sfg, log.p=TRUE)+dnorm((bg-v)/sbg, log=TRUE))
	sqrt(sfg^2+sbg^2) / sbg * integrate(normalConvolution,
        ifelse((bg-4*sbg)<0, 0, bg-4*sbg),
        bg+4*sbg, fg=fg, bg=bg, sfg=sfg,
        sbg=sbg, subdivisions=10000)$value
}
