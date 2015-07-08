read.idat <- function(idatfiles, bgxfile, dateinfo=FALSE, tolerance=0)
# Read idat data from gene expression BeadChips
# Matt Ritchie
# Created 30 September 2013. Modified on 14 January 2015 and 26 June 2015.
{
  if(!requireNamespace("illuminaio",quietly=TRUE)) stop("illuminaio package required but is not available")
  cat("Reading manifest file", bgxfile, "... ")
  bgx <- illuminaio::readBGX(bgxfile)
  cat("Done\n")
  nregprobes <-nrow(bgx$probes)
  nctrlprobes <-nrow(bgx$control)
  nprobes <- nregprobes+nctrlprobes
  nsamples <- length(idatfiles)
  elist <- new("EListRaw")
  elist$source <- "illumina"
  elist$genes <- rbind(bgx$probes[,c("Probe_Id","Array_Address_Id")], bgx$controls[,c("Probe_Id","Array_Address_Id")])
  elist$genes$Status <- "regular"
  elist$genes$Status[(nregprobes+1):nprobes] <- bgx$controls[,"Reporter_Group_Name"]
  elist$targets <- data.frame("IDATfile"=idatfiles, "DecodeInfo"=rep(NA, nsamples), "ScanInfo"=rep(NA, nsamples))
  tmp <- matrix(NA, nprobes, nsamples)
  colnames(tmp) <- idatfiles
  rownames(tmp) <- elist$genes[,"Array_Address_Id"]  
  elist$E <- elist$other$STDEV <- elist$other$NumBeads <- tmp
  rm(tmp)
  for(i in 1:nsamples) {
    cat("\t", idatfiles[i], "... ")
    tmp <- illuminaio::readIDAT(idatfiles[i])
    cat("Done\n")
    ind <- match(elist$genes[,"Array_Address_Id"], tmp$Quants[,'IllumicodeBinData'])
    if(sum(is.na(ind))>tolerance)
      stop("Can't match all ids in manifest with those in idat file ", idatfiles[i], "\n", sum(is.na(ind)), 
            " missing - please check that you have the right files, or consider setting \'tolerance\'=", sum(is.na(ind)))
    elist$E[!is.na(ind),i] <- round(tmp$Quants[ind[!is.na(ind)],'MeanBinData'],1) # intensity data
    elist$other$STDEV[!is.na(ind),i] <- tmp$Quants[ind[!is.na(ind)],'DevBinData'] # Bead STDEV
    elist$other$NumBeads[!is.na(ind),i] <- tmp$Quants[ind[!is.na(ind)],'NumGoodBeadsBinData'] # NumBeads
    if(dateinfo) {
      elist$targets[i,"DecodeInfo"] = paste(tmp$RunInfo[1,], sep="", collapse=" ")
      elist$targets[i,"ScanInfo"] = paste(tmp$RunInfo[2,], sep="", collapse=" ")
    }
  }
  cat("Finished reading data.\n")
  return(elist)
}
