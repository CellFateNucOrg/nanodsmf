


#' Split sequences with multiple motifs
#' @param tsv A tab serparated values text file with results from Nanopore CpG or GpC methylation calls
#' @param motif The motif to search for: "CG" or "GC".
#' @return A tsv with all mutli-motif sequences split into individual motifs with normalised log_lik_ratio
#' @examples
#' splitMotifs(MSssI_CpG,"CG")
#' @export
splitMotifs<-function(tsv,motif){
  newtsv<-tsv[rep(seq_len(nrow(tsv)), tsv$num_motifs),]
  motifmatch<-unlist(Biostrings::vmatchPattern(motif,tsv$sequence))
  newtsv$end<-newtsv$start+IRanges::end(motifmatch)
  newtsv$start<-newtsv$start+IRanges::start(motifmatch)
  newtsv$log_lik_ratio<-newtsv$log_lik_ratio/newtsv$num_motifs
  return(newtsv)
}


#' Convert nanopore tsv to methylation matrix
#' @param tsv A tab serparated values text file where individual motifs have been split
#' @param genomeGRs Genomic Ranges object for the regions to be analysed
#' @param motif Motif ("CG" or "GC" to for which the tsv was called)
#' @param binarise Convert log likelihoods to binary values: methylated(\eqn{ln(L) \ge 2.5}): 1; unmethylated(\eqn{ln(L) \le -2.5}): 0; inconclusive(\eqn{-2.5 < ln(L) < 2.5}): NA.  (default: binarise=TRUE)
#' @return A methylation matrix (reads x motif positions) with binary or log likelihood values
#' @examples
#' tsvToMethMat(splitMotifs(MSssI_CpG,"CG"),ttTi5605gr)
#' @export
tsvToMethMat<-function(tsv, genomeGRs, motif, binarise=TRUE){
  # make GR from tsv to subset only regions of interest
  tsvGR<-GenomicRanges::GRanges(seqnames=tsv$chromosome,
                                IRanges::IRanges(start=tsv$start,end=tsv$end),
                                strand=tsv$strand)
  keepCols<-c("chromosome","strand","start","read_name","log_lik_ratio" )
  GenomicRanges::mcols(tsvGR)<-tsv[,keepCols]
  matList=list()
  for (i in seq_along(genomeGRs)) {
    # deal with Cs on opposite strands
    ol<-GenomicRanges::findOverlaps(tsvGR,genomeGRs[i],ignore.strand=T)
    methGR<-tsvGR[S4Vectors::queryHits(ol)]
    options(tibble.width = Inf)
    methTab<-tidyr::spread(tibble::as.tibble(GenomicRanges::mcols(methGR)),
                           key=start,value=log_lik_ratio)
    methMat<-as.matrix(methTab[,-c(1,2,3)])
    rownames(methMat)<-methTab$read_name
    if (binarise==TRUE){
      methMat[abs(methMat) < 2.5]<- NA
      methMat[methMat >= 2.5]<- 1
      methMat[methMat <= -2.5]<- 0
    }
    matList[[genomeGRs[i]$ID]]<-methMat
  }
  return(matList)
}



