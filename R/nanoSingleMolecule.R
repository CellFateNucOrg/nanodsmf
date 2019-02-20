#' @importFrom magrittr "%>%"


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
#' @param genomeGR Genomic Ranges object for the regions to be analysed
#' @return A methylation matrix (reads x motif positions)
#' @examples
#' tsv2methMat(splitMotifs(MSssI_CpG,"CG"),genomeGR)
#' @export
tsv2methMat<-function(tsv,genomeGR,motif){
  # make GR from tsv to subset only regions of interest
  tsvGR<-GenomicRanges::GRanges(seqnames=tsv$chromosome,
                                IRanges::IRanges(start=tsv$start,end=tsv$end),
                                strand=tsv$strand)
  # deal with Cs on opposite strands
  keepCols<-c("chromosome","strand","start","read_name","log_lik_ratio" )
  GenomicRanges::mcols(tsvGR)<-tsv[,keepCols]
  ol<-GenomicRanges::findOverlaps(tsvGR,genomeGR)
  tsvGR<-tsvGR[S4Vectors::queryHits(ol)]
  options(tibble.width = Inf)
  tibble::as.tibble(GenomicRanges::mcols(tsvGR)) %>%
    dplyr::group_by(read_name) %>%
    tidyr::spread(key=start,value=log_lik_ratio)
  return(methMat)
}


