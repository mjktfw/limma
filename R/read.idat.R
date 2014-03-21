read.idat <- function(idatfiles, bgxfile, dateinfo=FALSE)
# Function to read idat data from gene expression BeadArrays
# Matt Ritchie, 30 September 2013
{
  require(illuminaio)
  cat("Reading manifest file", bgxfile, "... ")
  bgx <- readBGX(bgxfile)
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
    tmp <- readIDAT(idatfiles[i])
    cat("Done\n")
    ind <- match(elist$genes[,"Array_Address_Id"], tmp$Quants[,'IllumicodeBinData'])
    if(sum(is.na(ind))>0)
      stop("Can't match ids in manifest with those in idat file", idatfiles[i], "- please check that you have the right files\n")
    elist$E[,i] <- round(tmp$Quants[ind,'MeanBinData'],1) # intensity data
    elist$other$STDEV[,i] <- tmp$Quants[ind,'DevBinData'] # Bead STDEV
    elist$other$NumBeads[,i] <- tmp$Quants[ind,'NumGoodBeadsBinData'] # NumBeads
    if(dateinfo) {
      elist$targets[i,"DecodeInfo"] = paste(tmp$RunInfo[1,], sep="", collapse=" ")
      elist$targets[i,"ScanInfo"] = paste(tmp$RunInfo[2,], sep="", collapse=" ")
    }
  }
  cat("Finished reading data.\n")
  return(elist)
}
