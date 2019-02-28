


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
  df<-data.frame(read_name=newtsv$read_name,gnmstartfirst=newtsv$start,substrmatch=IRanges::start(motifmatch))
  #start in the tsv refers to position of first CG, not start of sequence that is in the sequence
  #column.
  df<-df %>%
    dplyr::group_by(read_name) %>%
    dplyr::mutate(minsubstr=min(substrmatch)) %>%
    dplyr::mutate(gnmstart=gnmstartfirst-minsubstr+substrmatch+1)
  # The first match is actuall always at position 6. but will keep this method just in case
  # nanopolish changes this
  newtsv$start<-df$gnmstart
  newtsv$end<-df$gnmstart+1
  newtsv$log_lik_ratio<-newtsv$log_lik_ratio/newtsv$num_motifs
  return(newtsv)
}


#' Convert nanopore tsv to genomic ranges
#' @param tsv A tab serparated values text file where individual motifs have been split
#' @param keepCols Columns from the tsv to put in as metadata in the new GR object
#' @return A genomic ranges object of motifs in tsv with some of the tsv data columns as metadata
#' @examples
#' tsvToGR(splitMotifs(MSssI_CpG,"CG"))
#' @export
tsvToGR<-function(tsv,keepCols=c("read_name","log_lik_ratio")){
  # make GR from tsv to subset only regions of interest
  tsvGR<-GenomicRanges::GRanges(seqnames=tsv$chromosome,
                                IRanges::IRanges(start=tsv$start,end=tsv$end),
                                strand=tsv$strand)
  GenomicRanges::mcols(tsvGR)<-tsv[,keepCols]
  return(tsvGR)
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
  tsvGR<-tsvToGR(tsv)
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



#' Merge CG tsv and GC tsv
#' @param tsvCG A tab serparated values text file where individual CG motifs have been split
#' @param tsvGC A tab serparated values text file where individual GC motifs have been split
#' @param genome A BSgenome, DNAstring or path to fasta file for the genome of interest
#' @return A tsv for both motifs with three contexts: HCG, GCH, GCGorCGC
#' @examples
#' mergeCGandGCtsv(splitMotifs(MSssI_CpG,"CG"),splitMotifs(MSssI_GpC,"GC"),ttTi5605dna)
#' @export
mergeCGandGCtsv<-function(tsvCG,tsvGC,genome){
  gnmMotifs<-findGenomeMotifs(genome)
  # make GR from tsvs
  tsvCGgr<-tsvToGR(tsvCG)
  tsvGCgr<-tsvToGR(tsvGC)
  # get isolated CGs
  ol<-GenomicRanges::findOverlaps(tsvCGgr,gnmMotifs[gnmMotifs$context=="HCG"])
  CGonly<-tsvCGgr[S4Vectors::queryHits(ol)]
  CGonly$context<-"HCG"
  # get isolated GCs
  ol<-GenomicRanges::findOverlaps(tsvGCgr,gnmMotifs[gnmMotifs$context=="GCH"])
  GConly<-tsvGCgr[S4Vectors::queryHits(ol)]
  GConly$context<-"GCH"
  # get CGs and GCs in runs
  ol<-GenomicRanges::findOverlaps(tsvCGgr,gnmMotifs[gnmMotifs$context=="GCGorCGC"])
  CGinRun<-tsvCGgr[S4Vectors::queryHits(ol)]
  ol<-GenomicRanges::findOverlaps(tsvGCgr,gnmMotifs[gnmMotifs$context=="GCGorCGC"])
  GCinRun<-tsvGCgr[S4Vectors::queryHits(ol)]
}




#' Merge CG and GC methylation matrices
#' @param CGmat A tab serparated values text file where individual CG motifs have been split
#' @param GCmat A tab serparated values text file where individual GC motifs have been split
#' @param genome A BSgenome, DNAstring or path to fasta file for the genome of interest
#' @return A tsv for both motifs with three contexts: HCG, GCH, GCGorCGC
#' @examples
#' mergeCGandGCtsv(splitMotifs(MSssI_CpG,"CG"),splitMotifs(MSssI_GpC,"GC"),ttTi5605dna)
#' @export
mergeCGandGCtsv<-function(tsvCG,tsvGC,genome){
  gnmMotifs<-findGenomeMotifs(genome)
  # make GR from tsvs
  tsvCGgr<-tsvToGR(tsvCG)
  tsvGCgr<-tsvToGR(tsvGC)
  # get isolated CGs
  ol<-GenomicRanges::findOverlaps(tsvCGgr,gnmMotifs[gnmMotifs$context=="HCG"])
  CGonly<-tsvCGgr[S4Vectors::queryHits(ol)]
  CGonly$context<-"HCG"
  # get isolated GCs
  ol<-GenomicRanges::findOverlaps(tsvGCgr,gnmMotifs[gnmMotifs$context=="GCH"])
  GConly<-tsvGCgr[S4Vectors::queryHits(ol)]
  GConly$context<-"GCH"
  # get CGs and GCs in runs
  ol<-GenomicRanges::findOverlaps(tsvCGgr,gnmMotifs[gnmMotifs$context=="GCGorCGC"])
  CGinRun<-tsvCGgr[S4Vectors::queryHits(ol)]
  ol<-GenomicRanges::findOverlaps(tsvGCgr,gnmMotifs[gnmMotifs$context=="GCGorCGC"])
  GCinRun<-tsvGCgr[S4Vectors::queryHits(ol)]
  runs<-c(CGinRun,GCinRun)
  for (readName in unique(runs$read_name))
    applyGRonGR(gnmMotifs,runs[runs$read_name==readName],"log_lik_ratio",sum)
}

