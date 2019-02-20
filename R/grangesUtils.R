# general useful functions for granges

#' Convert DNAStringSet to GenomicRanges
#' @param DNAss A DNAStringSet
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' DNAstringset2gr(Biostrings::DNAStringSet(x=c(chr1="agctgtagct",chr2="agagagagttt"))
#' @export
DNAstringset2gr<-function(DNAss,defaultStrand="*") {
  #extract only the first word in fasta '>' row to be the seq name
  chrNames<-sapply(strsplit(names(DNAss)," "),'[[',1)
  gr<-GenomicRanges::GRanges(seqnames=chrNames,
                             ranges=IRanges::IRanges(start=1,width=Biostrings::width(DNAss)),
                             strand=defaultStrand)
  return(gr)
}


#' importFrom magrittr "%>%"

#' Convert DNAStringSet to GenomicRanges
#' @param BSgenome A BSgenome object
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' BSgenome2gr(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenome2gr<-function(BSgnm,defaultStrand="*") {
  gr<-GenomicRanges::GRanges(seqnames=BSgenome::seqnames(BSgnm),
                             ranges=IRanges::IRanges(start=1,end=GenomeInfoDb::seqlengths(BSgnm)),
                             strand="*")
  return(gr)
}



#' Convert UCSC to wormbase chromosome names
#' @param ucscGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' ucsc2wbGR(BSgenome2gr(BSgenome.Celegans.UCSC.ce11::Celegans))
#' @export
ucsc2wbGR<-function(ucscGR) {
  wbGR<-ucscGR
  GenomeInfoDb::seqlevels(wbGR)<-gsub("chr","",GenomeInfoDb::seqlevels(wbGR))
  GenomeInfoDb::seqlevels(wbGR)<-gsub("M","MtDNA",GenomeInfoDb::seqlevels(wbGR))
  return(wbGR)
}



#' Convert wormbase to ucsc chromosome names
#' @param wbGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' wb2ucscGR(ucsc2wbGR(BSgenome2gr(BSgenome.Celegans.UCSC.ce11::Celegans)))
#' @export
wb2ucscGR<-function(wbGR) {
  ucscGR<-wbGR
  GenomeInfoDb::seqlevels(ucscGR)<-gsub("MtDNA","M",GenomeInfoDb::seqlevels(ucscGR))
  GenomeInfoDb::seqlevels(ucscGR)<-past0("chr",GenomeInfoDb::seqlevels(ucscGR))
  return(ucscGR)
}

