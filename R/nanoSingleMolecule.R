


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
#' @param motif Motif ("CG" or "GC" to for which the tsv was called)
#' @param binarise Convert log likelihoods to binary values: methylated(\eqn{ln(L) \ge 2.5}): 1; unmethylated(\eqn{ln(L) \le -2.5}): 0; inconclusive(\eqn{-2.5 < ln(L) < 2.5}): NA.  (default: binarise=TRUE)
#' @return A methylation matrix (reads x motif positions) with binary or log likelihood values
#' @examples
#' tsvToMethMat(splitMotifs(MSssI_CpG,"CG"),genomeGR)
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





#' Do a single molecule plot of dSMF data
#' @param mat A methylation matrix
#' @param regionName The name of the region to be plotted (usually gene name). Must match Ids in reagionGRs
#' @param regionGRs Genomic Ranges which include region to be plotted (ID metadata column must match regionName)
#' @param featureGRs Genomic Ranges of any features of interest (e.g. TSS) that you wish to be plotted as vertical red line
#' @param myXlab Text to be used for X axis label (default="CpG/GpC position")
#' @param featureLabel  Text to be used to label special feature you are plotting (default="TSS")
#' @param title Text to be used as plot title (default=NULL will results in the plot named accorind to the regionGR)
#' @param baseFontSize The base font size to be used for the plot (default is 12)
#' @return A ggplot2 plotting object
#' @examples
#' plotSingleMoleculesAmp(methMat,"ttTi5605",genomeGR,featureGRs=c(),myXlab="CpG position",featureLabel="TSS",title=NULL,baseFontSize=12)
#' @export
plotSingleMoleculesAmp<-function(mat,regionName,regionGRs,featureGRs=c(),myXlab="CpG/GpC position",featureLabel="TSS",title=NULL,
                                 baseFontSize=12) {
  ### single molecule plot. mat is matrix containing methylation values at different postions (columns) in
  # individual reads (rows). regionName is the ID of the amplicon or genomic region being plotted. regionGRs is a
  # genomicRanges object containing the region being plotted. one of its mcols must have a name "ID" in which the
  # same ID as in regionName appears. featureGRs is genomic ranges object for plotting location of some feature in
  # the region, such as the TSS. myXlab is the X axis label. featureLabel is the label for the type of feature that
  # will be plotted underneath the feature
  if(dim(mat)[1]>10) {
    regionGR<-regionGRs[match(regionName,regionGRs$ID)]
    if (length(featureGRs)>0) {
      featGR<-featureGRs[match(regionName,featureGRs$ID)]
    }
    na.matrix<-is.na(mat)
    mat[na.matrix]<- -1
    # try to perform heirarchical clustering
    hc <- try(
      hclust(dist(apply(mat,2,as.numeric))),
      silent = TRUE)
    mat[na.matrix]<-NA
    if (class(hc) == "try-error") {
      df<-as.data.frame(mat)
      print("hclust failed. Matrix dim: ")
      print(dim(mat))
    } else {
      df<-as.data.frame(mat[hc$order,])
    }

    reads<-row.names(df)
    d<-tidyr::gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      title=paste0(regionName, ": ", GenomicRanges::seqnames(regionGR)," ", GenomicRanges::strand(regionGR),"ve strand")
    }
    p<-ggplot(d,aes(x=position,y=molecules,width=2)) +
      geom_tile(aes(width=6*GenomicRanges::width(regionGR)/500,fill=methylation),alpha=0.8) +
      scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white",
                        labels=c("protected","accessible"),name="dSMF") + theme_light(base_size=baseFontSize) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold",hjust = 0.5),
            legend.position="bottom", legend.key.height = unit(0.1, "cm")) +
      ggtitle(title) +
      xlab(myXlab) + ylab("Single molecules") + xlim(start(regionGR),end(regionGR)+10)
    if(length(featureGRs)>0) {
      p<-p+geom_linerange(aes(x=start(featGR), y=NULL, ymin=0, ymax=length(reads)+max(3,0.04*length(reads))),col="red") +
        ggplot2::annotate("segment", x = GenomicRanges::start(featGR), xend = GenomicRanges::start(featGR)+20*ifelse(GenomicRanges::strand(featGR)=="-",-1,1),
                          y = length(reads)+max(3,0.04*length(reads)), yend =length(reads)+max(3,0.04*length(reads)), colour = "red",
                          arrow=arrow(length = unit(0.3, "cm")), size=0.7) +
        ggplot2::annotate(geom="text", x=GenomicRanges::start(featGR), y=-max(2,0.03*length(reads)), label=featureLabel,color="red")
    }
  } else {
    p<-NULL
  }
  return(p)
}



#' Find non overlapping CG, GC and GCG/CGC motifs in genome
#' @param genomeFile A DNAstringSet for the genome
findGenomeMotifs<-function(genome){

  # deal with different input genome formats
  if (class(genome)=="BSgenome"){
    genome<-BSgenomeToDNAStringSet(genome)
  } else if (is.character(genome)){
    genome<-readDNAStringSet(genome)

  }
  #strip of anything after a space to deal with additional fasta header info
  names(genome)<-gsub("[[:space:]].*$","",names(genome),perl=F)

  gnmCGs<-mIdxToGR(Biostrings::vmatchPattern("CG",genome))
  gnmGCs<-mIdxToGR(Biostrings::vmatchPattern("GC",genome))
  gnmGCGs<-mIdxToGR(Biostrings::vmatchPattern("GCG",genome))
  gnmCGCs<-mIdxToGR(Biostrings::vmatchPattern("CGC",genome))




}

findNonOverlappingMotifs<-function(DNAss) {
  gnmCGs<-Biostrings::vmatchPattern("CG",genome)
  gnmGCs<-Biostrings::vmatchPattern("GC",genome)
  gnmGCGs<-Biostrings::vmatchPattern("GCG",genome)
  gnmCGCs<-Biostrings::vmatchPattern("CGC",genome)
}


